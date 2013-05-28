#include "wiggleIterators.h"
#include "sam.h"

typedef struct bamReaderData_st {
	mplp_aux_t *data;
	int n_plp, tid0 = -1, ref_len, ref_tid = -1, max_depth, max_indel_depth;
	const bam_pileup1_t *plp;
	bam_mplp_t iter;
	bam_header_t *h = 0;
	void *rghash = 0;
	bam_sample_t *sm = 0;
	kstring_t buf;
	mplp_pileup_t gplp;
} BamReaderData;

void OpenBamFile(BamFileReader * data, char * filename) {
	bam_header_t *h_tmp;

	memset(&gplp, 0, sizeof(mplp_pileup_t));
	memset(&buf, 0, sizeof(kstring_t));
	memset(&bc, 0, sizeof(bcf_call_t));
	data = calloc(n, sizeof(void*));
	plp = calloc(n, sizeof(void*));
	n_plp = calloc(n, sizeof(int*));
	sm = bam_smpl_init();
	data = calloc(1, sizeof(mplp_aux_t));

	// read the header and initialize data
	if (strcmp(fiename, "-"))
		data->fp = bam_open(filename, "r");
	else
		data->fp = bam_dopen(fileno(stdin), "r");
	data->conf = conf;
	h_tmp = bam_header_read(data->fp);
	data->h = h_tmp; 
	bam_smpl_add(sm, filename, (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
	rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);
	int beg, end;
	bam_index_t *idx;
	idx = bam_index_load(filename);
	if (idx == 0) {
		fprintf(stderr, "[%s] fail to load index for input.\n", __func__);
		exit(1);
	}
	if (bam_parse_region(h_tmp, conf->reg, &tid, &beg, &end) < 0) {
		fprintf(stderr, "[%s] malformatted region or wrong seqname for input.\n", __func__);
		exit(1);
	}
	tid0 = tid, beg0 = beg, end0 = end;
	data->iter = bam_iter_query(idx, tid, beg, end);
	bam_index_destroy(idx);

	gplp.n = sm->n;
	gplp.n_plp = calloc(sm->n, sizeof(int));
	gplp.m_plp = calloc(sm->n, sizeof(int));
	gplp.plp = calloc(sm->n, sizeof(void*));

	fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
	if (tid0 >= 0 && conf->fai) { // region is set
		ref = faidx_fetch_seq(conf->fai, h_tmp->target_name[tid0], 0, 0x7fffffff, &ref_len);
		ref_tid = tid0;
		data->ref = ref;
		data->ref_id = tid0;
	} else ref_tid = -1, ref = 0;
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	max_depth = conf->max_depth;
	if (max_depth * sm->n > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
	if (max_depth * sm->n < 8000) {
		max_depth = 8000 / sm->n;
		fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
	}
	max_indel_depth = conf->max_indel_depth * sm->n;
	bam_mplp_set_maxcnt(iter, max_depth);

}

void CloseBamFile(BamReaderData * data) {
	// Seriously, does samtools not provide any convience destructors!??
	bam_smpl_destroy(data->sm); 
	free(data->buf.s);
	free(data->gplp.plp); 
	free(data->gplp.n_plp); 
	free(data->gplp.m_plp);
	bam_mplp_destroy(data->iter);
	bam_header_destroy(data->h);
	bam_close(data->data->fp);
	if (data->data->iter) 
		bam_iter_destroy(data->data->iter);
	free(data->data); 
	free(data->plp); 
	free(data->n_plp);
}

void BamReaderPop(WiggleIterator * wi)
{
	BamReaderData * data = (BamReaderData *) wi->data;
	int j, tid, n_pos, cnt;

	if (bam_mplp_auto(data->iter, &tid, &pos, &data->n_plp, &data->plp) > 0) {
		// If fell out of bounds:
		if (wi->stop > 0 && pos >= wi->stop) {
			wi->done = true;
			return;
		} 
		
		// If a new chrom. was entered:
		if (tid != data->ref_tid) {
			if (data->conf->fai) 
				data->nextChrom = faidx_fetch_seq(data->conf->fai, data->h->target_name[tid], 0, 0x7fffffff, &data->ref_len);
			else
				data->nextChrom = NULL;
			data->data->ref = data->nextChrom;
			data->data->ref_id = tid;
			data->ref_tid = tid;
		}

		// Count reads in pileup
		cnt = 0;
		const bam_pileup1_t *p = data->plp;
		for (j = 0; j < n_plp; ++j) {
			if (bam1_qual(p->b)[p->qpos] >= data->conf->min_baseQ) 
				cnt++;
			p++;
		}

		// Its a wrap:
		// +1 to account for 0-based indexing in BAMs:
		wi->nextStart = pos + 1;
		wi->nextFinish = pos + 2;
		wi->nextValue = cnt;
	} else {
		closeBamFile(data);
		wi->done = true;
	}
}

void BamReaderSeek(WiggleIterator * wi, char * chrom, int start, int finish) {
	BamReaderData * data = (BamReaderData *) wi->data;
	wi->chrom = chrom;
	// TODO
	wi->stop = finish;
}

WiggleIterator * BamReader(char * filename) {
	BamReaderData * data = (BamReaderData *) calloc(1, sizeof(BamReaderData));
	OpenBamFile(data, filename);
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek);
}
