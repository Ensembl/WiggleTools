// Copyright [1999-2017] EMBL-European Bioinformatics Institute
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

void libBigWigInit(int size) {
	if(bwInit(size) != 0) {
		fprintf(stderr, "Received an error in bwInit\n");
		exit(1);
	}
}

static int readIteratorIntervals(bwOverlapIterator_t *iter, char * chrom, int stretch_start, int stretch_stop, BigWiggleReaderData * data) {
	int index;
	for(index = 0; index < iter->intervals->l; index++) {
		int start = iter->intervals->start[index] + 1;
		int finish = iter->intervals->end[index] + 1;

		// Box into queried stretch
		start = start < stretch_start? stretch_start: start;
		finish = finish < stretch_stop? finish: stretch_stop;

		if (pushValuesToBuffer(data->bufferedReaderData, chrom, start, finish, iter->intervals->value[index]))
			return 1;
	}
	return 0;
}

static int readBigWiggleRegion(BigWiggleReaderData * data, char * chrom, int start, int stop) {
	// This hack is required because libBigWig does not handle negative numbers
	if (start < 1)
		start = 1;
	if (stop < 1)
		stop = 1;
	bwOverlapIterator_t *iter = bwOverlappingIntervalsIterator(data->fp, chrom, start - 1, stop - 1, MAX_BLOCKS);
	if (!iter)
		return 0;

	while(iter->data) {
		if (readIteratorIntervals(iter, chrom, start, stop, data)) {
			bwIteratorDestroy(iter);
			return 1;
		}
		iter = bwIteratorNext(iter);
	}
	bwIteratorDestroy(iter);
	return 0;
}

static int readBigWiggleChromosome(BigWiggleReaderData * data, char * chrom, int length) {
	int start;
	int stretch=10000;

	for (start = 1; start < length; start+=stretch) {
		if (readBigWiggleRegion(data, chrom, start, start+stretch))
			return 1;
	}

	return 0;
}

void * readBigWiggle(void * ptr) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) ptr;
	if (data->chrom)
		readBigWiggleRegion(data, data->chrom, data->start, data->stop);
	else {
		int chrom_index;
		Chrom_length * chrom_lengths = calloc(data->fp->cl->nKeys, sizeof(Chrom_length));
		for (chrom_index = 0; chrom_index < data->fp->cl->nKeys; chrom_index++) {
			chrom_lengths[chrom_index].chrom = data->fp->cl->chrom[chrom_index];
			chrom_lengths[chrom_index].length = data->fp->cl->len[chrom_index];
		}

		qsort(chrom_lengths, data->fp->cl->nKeys, sizeof(Chrom_length), compare_chrom_lengths);

		for (chrom_index = 0; chrom_index < data->fp->cl->nKeys; chrom_index++)
			if (readBigWiggleChromosome(data, chrom_lengths[chrom_index].chrom, chrom_lengths[chrom_index].length))
				break;

		free(chrom_lengths);
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
		printf("File %s is not in BigWig format\n", filename);
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
	return newWiggleIterator(data, &BigWiggleReaderPop, &BigWiggleReaderSeek, 0, false);
}	
