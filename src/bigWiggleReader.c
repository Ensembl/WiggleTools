// Copyright [1999-2016] EMBL-European Bioinformatics Institute
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
#include "bigWig.h"
#include "bufferedReader.h"

static int MAX_BLOCKS = 100;

typedef struct bigWiggleReaderData_st {
	bigWigFile_t * fp;
	char * chrom;
	int start;
	int stop;
	BufferedReaderData * bufferedReaderData;
} BigWiggleReaderData;

static int readIteratorIntervals(bwOverlapIterator_t *iter, char * chrom, BigWiggleReaderData * data) {
	int index;
	for(index = 0; index < iter->intervals->l; index++) {
		int start = iter->intervals->start[index] + 1;
		int finish = iter->intervals->end[index] + 1;
		if (data->stop > 0) {
			if (start >= data->stop)
				return 1;
			else if (finish > data->stop)
				finish = data->stop;
		}
		if (pushValuesToBuffer(data->bufferedReaderData, chrom, start, finish, iter->intervals->value[index]))
			return 1;
	}
	return 0;
}

static int readBigWiggleRegion(BigWiggleReaderData * data, char * chrom, int start, int stop) {
	bwOverlapIterator_t *iter = bwOverlappingIntervalsIterator(data->fp, chrom, start, stop, MAX_BLOCKS);
	while(iter->data) {
		if (readIteratorIntervals(iter, chrom, data))
			return 1;
		iter = bwIteratorNext(iter);
	}
	bwIteratorDestroy(iter);
	return 0;
}

void * readBigWiggle(void * ptr) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) ptr;
	if (data->chrom)
		readBigWiggleRegion(data, data->chrom, data->start, data->stop);
	else {
		int chrom_index;
		for (chrom_index = 0; chrom_index < data->fp->cl->nKeys; chrom_index++)
			if (readBigWiggleRegion(data, data->fp->cl->chrom[chrom_index], 0, data->fp->cl->len[chrom_index]))
				break;
	}

	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void BigWiggleReaderPop(WiggleIterator * wi) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void openBigWiggle(BigWiggleReaderData * data, char * filename, bool holdFire) {
	if(!bwIsBigWig(filename, NULL)) {
		printf("File %s is not in BigBEd format", filename);
		exit(1);
	}
	data->fp = bwOpen(filename, NULL, "r");
	if (!holdFire)
		launchBufferedReader(&readBigWiggle, data, &(data->bufferedReaderData));
}

void BigWiggleReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) wi->data; 

	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}
	data->chrom = chrom;
	data->start = start;
	data->stop = finish;
	launchBufferedReader(&readBigWiggle, data, &(data->bufferedReaderData));
	wi->done = false;
	BigWiggleReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigWiggleReaderPop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
}

WiggleIterator * BigWiggleReader(char * f, bool holdFire) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) calloc(1, sizeof(BigWiggleReaderData));
	openBigWiggle(data, f, holdFire);
	return newWiggleIterator(data, &BigWiggleReaderPop, &BigWiggleReaderSeek, 0);
}	
