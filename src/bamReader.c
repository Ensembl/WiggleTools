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
#include <pthread.h>

#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80

static int MAX_HEAD_START = 30;
static int BLOCK_SIZE = 10000;
extern int mplp_func(void *data, bam1_t *b);

typedef struct blockData_st {
	char **chrom;
	int * start;
	int * value;
	int index;
	int count;
	struct blockData_st * next;
} BlockData;

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
	pthread_t downloaderThreadID;

	// Output of downloader
	BlockData * blockData, *lastBlockData;
	pthread_mutex_t count_mutex;
	pthread_cond_t count_cond;
	int blockCount;

	// BAM stuff
	mplp_conf_t * conf;
	bam_index_t * idx;
	mplp_aux_t * data;
	int ref_tid;
	bam_mplp_t iter;
} BamReaderData;

static BlockData * createBlockData() {
	BlockData * new = (BlockData * ) calloc(1, sizeof(BlockData));
	new->chrom = (char **) calloc(BLOCK_SIZE, sizeof(char*));
	new->start = (int *) calloc(BLOCK_SIZE, sizeof(int));
	new->value = (int *) calloc(BLOCK_SIZE, sizeof(int));
	return new;
}

void destroyBamBlockData(BlockData * data) {
	free(data->chrom);
	free(data->start);
	free(data->value);
	free(data);
}

static void waitForNextBlock(BamReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	// Check whether allowed to step forward
	if (data->blockCount == 0) {
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	}
	// Signal freed memory
	data->blockCount--;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
}

void goToNextBamBlock(BamReaderData * data) {
	BlockData * prevBlockData = data->blockData;
	data->blockData = data->blockData->next;
	destroyBamBlockData(prevBlockData);
}

static bool declareNewBlock(BamReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);

	if (data->blockCount > MAX_HEAD_START) {
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	} 
	if (data->blockCount < 0) {
		pthread_mutex_unlock(&data->count_mutex);
		return true;
	}
	data->blockCount++;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
	return false;
}

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

		if (data->blockData == NULL) {
			data->lastBlockData = data->blockData = createBlockData();
			//declareNewBlock(data);
		} else if (data->lastBlockData->count == BLOCK_SIZE) {
			data->lastBlockData->next = createBlockData();
			data->lastBlockData = data->lastBlockData->next;
			declareNewBlock(data);
		}

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

		if (data->stop > 0 && (strcmp(chrom, data->chrom) == 0 && pos >= data->stop)) {
			break;
		}

		int index = data->lastBlockData->count;
		data->lastBlockData->chrom[index] = chrom;
		// +1 to account for 0-based indexing in BAMs:
		data->lastBlockData->start[index] = pos + 1;
		data->lastBlockData->value[index] = cnt;
		data->lastBlockData->count++;
	}

	// Signals to the reader that it can step into NULL block, hence marking the end of the download
	declareNewBlock(data);
	declareNewBlock(data);
	return NULL;
}

void launchBamDownloader(BamReaderData * data) {
	pthread_mutex_init(&data->count_mutex, NULL);
	pthread_cond_init(&data->count_cond, NULL);

	int err = pthread_create(&data->downloaderThreadID, NULL, &downloadBamFile, data);
	if (err) {
		printf("Could not create new thread %i\n", err);
		abort();
	}

	waitForNextBlock(data);
}

void killBamDownloader(BamReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	data->blockCount = -1;
	// Send a signal in case the slave is waiting somewhere
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
	pthread_join(data->downloaderThreadID, NULL);

	pthread_mutex_destroy(&data->count_mutex);
	pthread_cond_destroy(&data->count_cond);

	if (data->data->iter) 
		bam_iter_destroy(data->data->iter);

	while (data->blockData) {
		BlockData * prevData = data->blockData;
		data->blockData = data->blockData->next;
		destroyBamBlockData(prevData);
	}

	data->lastBlockData = NULL;
	data->blockCount = 0;
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

	if (wi->done)
		return;
	else if (data->blockData == NULL) {
		wi->done = true;
		return;
	} else if (data->blockData->index == data->blockData->count) {
		waitForNextBlock(data);
		goToNextBamBlock(data);
		if (data->blockData == NULL) {
			wi->done = true;
			return;
		}
	} 

	int index = data->blockData->index;
	wi->chrom = data->blockData->chrom[index];
	wi->start = data->blockData->start[index];
	wi->finish = data->blockData->start[index] + 1;
	wi->value = (double) data->blockData->value[index];
	data->blockData->index++;
}

void BamReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	char region[1000];

	killBamDownloader(wi->data);

	sprintf(region, "%s:%i-%i", chrom, start, finish);
	BamReaderData * data = (BamReaderData *) wi->data;
	if (data->conf->reg)
		free(data->conf->reg);
	data->conf->reg = region;
	seekRegion(data);
	launchBamDownloader(data);
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
	launchBamDownloader(data);
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek);
}
