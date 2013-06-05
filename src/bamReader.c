// Copyright (c) 2013, Daniel Zerbino
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
// (1) Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer. 
// 
// (2) Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in
// the documentation and/or other materials provided with the
// distribution.  
// 
// (3)The name of the author may not be used to
// endorse or promote products derived from this software without
// specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <string.h>
#include "sam.h"
#include "faidx.h"
#include "wiggleIterators.h"

#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80

extern int mplp_func(void *data, bam1_t *b);

typedef struct {
	int max_mq, min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
	int openQ, extQ, tandemQ, min_support; // for indels
	double min_frac; // for indels
	char *reg, *pl_list;
	faidx_t *fai;
	void *bed, *rghash;
} mplp_conf_t;

typedef struct {
	bamFile fp;
	bam_iter_t iter;
	bam_header_t *h;
	int ref_id;
	char *ref;
	const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
	int n;
	int *n_plp, *m_plp;
	bam_pileup1_t **plp;
} mplp_pileup_t;

typedef struct bamReaderData_st {
	mplp_conf_t * conf;
	bam_index_t * idx;
	mplp_aux_t * data;
	int ref_tid;
	bam_mplp_t iter;
} BamReaderData;

void setSamtoolsDefaultConf(BamReaderData * data) {
	data->conf = (mplp_conf_t *) calloc(1, sizeof(mplp_conf_t));
	memset(data->conf, 0, sizeof(mplp_conf_t));
	data->conf->max_mq = 60;
	data->conf->min_baseQ = 13;
	data->conf->capQ_thres = 0;
	data->conf->max_depth = 250; 
	data->conf->max_indel_depth = 250;
	data->conf->openQ = 40; 
	data->conf->extQ = 20; 
	data->conf->tandemQ = 100;
	data->conf->min_frac = 0.002; 
	data->conf->min_support = 1;
	data->conf->flag = MPLP_NO_ORPHAN | MPLP_REALN;
}

void seekRegion(BamReaderData * data) {
	int max_depth, tid, beg, end;

	if (data->conf->reg) {
		// Create BAM iterator at region
		if (bam_parse_region(data->data->h, data->conf->reg, &tid, &beg, &end) < 0) {
			fprintf(stderr, "[%s] malformatted region or wrong seqname for input.\n", __func__);
			exit(1);
		}
		data->data->iter = bam_iter_query(data->idx, tid, beg, end);
		data->ref_tid = tid;
	} else {
		// Create general BAM iterator
		data->data->iter = NULL;
	}

	// Create pileup iterator	
	data->iter = bam_mplp_init(1, mplp_func, (void**) &data->data);

	// Set max depth
	max_depth = data->conf->max_depth;
	if (max_depth > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
	if (max_depth < 8000) {
		max_depth = 8000;
		fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
	}
	bam_mplp_set_maxcnt(data->iter, max_depth);

}

void OpenBamFile(BamReaderData * data, char * filename) {
	// Allocate space
	data->data = (mplp_aux_t *) calloc(1, sizeof(mplp_aux_t));

	// read the header and initialize data
	if (strcmp(filename, "-"))
		data->data->fp = bam_open(filename, "r");
	else
		data->data->fp = bam_dopen(fileno(stdin), "r");
	data->data->conf = data->conf;
	data->data->h = bam_header_read(data->data->fp);

	// Load index
	data->idx = bam_index_load(filename);
	if (data->idx == 0) {
		fprintf(stderr, "[%s] fail to load index for input.\n", __func__);
		exit(1);
	}

	// Start reading
	data->ref_tid = -1;
	seekRegion(data);
}

void closeBamFile(BamReaderData * data) {
	// Seriously, does samtools not provide any convience destructors!??
	bam_mplp_destroy(data->iter);
	//bam_header_destroy(data->data->h);
	bam_close(data->data->fp);
	if (data->data->iter) 
		bam_iter_destroy(data->data->iter);
	free(data->data); 
	bam_index_destroy(data->idx);
}

void BamReaderPop(WiggleIterator * wi)
{
	BamReaderData * data = (BamReaderData *) wi->data;
	int j, tid, cnt, pos, n_plp;
	const bam_pileup1_t *plp;

	if (wi->done)
		return;

	if (bam_mplp_auto(data->iter, &tid, &pos, &n_plp, &plp) > 0) {
		// TODO Check you are in bounds
		// If a new chrom. was entered:
		if (tid != data->ref_tid) {
			wi->chrom = data->data->h->target_name[tid];
			data->ref_tid = tid;
		}

		// Count reads in pileup
		cnt = 0;
		const bam_pileup1_t *p = plp;
		for (j = 0; j < n_plp; ++j) {
			if (bam1_qual(p->b)[p->qpos] >= data->conf->min_baseQ) 
				cnt++;
			p++;
		}

		// Its a wrap:
		// +1 to account for 0-based indexing in BAMs:
		wi->start = pos + 1;
		wi->finish = pos + 2;
		wi->value = cnt;
	} else {
		closeBamFile(data);
		wi->done = true;
	}
}

void BamReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	char region[1000];
	sprintf(region, "%s:%i-%i", chrom, start, finish);
	BamReaderData * data = (BamReaderData *) wi->data;
	if (data->conf->reg)
		free(data->conf->reg);
	data->conf->reg = region;
	seekRegion(data);
}

WiggleIterator * BamReader(char * filename) {
	BamReaderData * data = (BamReaderData *) calloc(1, sizeof(BamReaderData));
	setSamtoolsDefaultConf(data);
	OpenBamFile(data, filename);
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek);
}
