// Copyright 2013 EMBL-EBI
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string.h>
#include "sam.h"
#include "faidx.h"
#include "wiggleIterator.h"
#include "bufferedReader.h"

#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80

extern int mplp_func(void *data, bam1_t *b);

typedef struct {
	int max_mq, min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter;
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


typedef struct bamFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int start, stop;
	BufferedReaderData * bufferedReaderData;

	// BAM stuff
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
	data->conf->min_baseQ = 0;
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

static void * downloadBamFile(void * args) {
	BamReaderData * data = (BamReaderData *) args;
	int j, tid, cnt, pos, n_plp;
	const bam_pileup1_t *plp;

	while (bam_mplp_auto(data->iter, &tid, &pos, &n_plp, &plp) > 0) {
		// Count reads in pileup
		cnt = 0;
		const bam_pileup1_t *p = plp;
		for (j = 0; j < n_plp; ++j) {
			//if (bam1_qual(p->b)[p->qpos] >= data->conf->min_baseQ) 
				cnt++;
			p++;
		}

		// Its a wrap:
		char * chrom = data->data->h->target_name[tid];

		if (data->stop > 0 && (strcmp(chrom, data->chrom) == 0 && pos >= data->stop))
			break;

		// +1 to account for 0-based indexing in BAMs:
		if (pushValuesToBuffer(data->bufferedReaderData, chrom, pos+1, pos+2, cnt))
			break;
	}

	bam_iter_destroy(data->data->iter);
	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void seekRegion(BamReaderData * data) {
	int tid, beg, end;

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

void BamReaderPop(WiggleIterator * wi) {
	BamReaderData * data = (BamReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void BamReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	char region[1000];

	killBufferedReader(wi->data);

	sprintf(region, "%s:%i-%i", chrom, start, finish);
	BamReaderData * data = (BamReaderData *) wi->data;
	if (data->conf->reg)
		free(data->conf->reg);
	data->conf->reg = region;
	seekRegion(data);
	launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	wi->done = false;
	BamReaderPop(wi);

	while (strcmp(wi->chrom, chrom) < 0 || wi->finish <= start)
		BamReaderPop(wi);

	data->chrom = chrom;
	data->stop = finish;
}

WiggleIterator * BamReader(char * filename) {
	BamReaderData * data = (BamReaderData *) calloc(1, sizeof(BamReaderData));
	setSamtoolsDefaultConf(data);
	OpenBamFile(data, filename);
	launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek);
}
