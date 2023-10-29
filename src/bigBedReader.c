// Copyright [1999-2017] EMBL-European Bioinformatics Institute

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

typedef struct bigBedReaderData_st {
	bigWigFile_t * fp;
	char * chrom;
	int start;
	int stop;
	BufferedReaderData * bufferedReaderData;
} BigBedReaderData;

static int readIteratorEntries(bwOverlapIterator_t *iter, char * chrom, int stretch_start, int stretch_stop, BigBedReaderData * data) {
	int index;
	for(index = 0; index < iter->entries->l; index++) {
		int start = iter->entries->start[index] + 1;
		int finish = iter->entries->end[index] + 1;

		// Box into queried stretch
		start = start < stretch_start? stretch_start: start;
		finish = finish < stretch_stop? finish: stretch_stop;

		if (pushValuesToBuffer(data->bufferedReaderData, chrom, start, finish, 1))
			return 1;
	}
	return 0;
}

static int readBigBedRegion(BigBedReaderData * data, char * chrom, int start, int stop) {
	// This hack is required because libBigWig does not handle negative numbers
	if (start < 1)
		start = 1;
	if (stop < 1)
		stop = 1;
	// BigBed format 1 indexed, hence the -1s
	bwOverlapIterator_t *iter = bbOverlappingEntriesIterator(data->fp, chrom, start - 1, stop - 1, 0, MAX_BLOCKS);
	if (!iter)
		return 0;
	
	while(iter->data) {
		if (readIteratorEntries(iter, chrom, start, stop, data)) {
			bwIteratorDestroy(iter);
			return 1;
		}
		iter = bwIteratorNext(iter);
	}
	bwIteratorDestroy(iter);
	return 0;
}

void * readBigBed(void * ptr) {
	BigBedReaderData * data = (BigBedReaderData *) ptr;

	if (data->chrom)
		readBigBedRegion(data, data->chrom, data->start, data->stop);
	else {
		int chrom_index;
		Chrom_length * chrom_lengths = calloc(data->fp->cl->nKeys, sizeof(Chrom_length));
		for (chrom_index = 0; chrom_index < data->fp->cl->nKeys; chrom_index++) {
			chrom_lengths[chrom_index].chrom = data->fp->cl->chrom[chrom_index];
			chrom_lengths[chrom_index].length = data->fp->cl->len[chrom_index];
		}

		qsort(chrom_lengths, data->fp->cl->nKeys, sizeof(Chrom_length), compare_chrom_lengths);

		for (chrom_index = 0; chrom_index < data->fp->cl->nKeys; chrom_index++)
			if (readBigBedRegion(data, chrom_lengths[chrom_index].chrom, 1, chrom_lengths[chrom_index].length))
				break;

		free(chrom_lengths);
	}

	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void BigBedReaderPop(WiggleIterator * wi) {
	BigBedReaderData * data = (BigBedReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void openBigBed(BigBedReaderData * data, char * filename, bool holdFire) {
	if(!bbIsBigBed(filename, NULL)) {
		printf("File %s is not in BigBed format\n", filename);
		exit(1);
	}
	data->fp = bbOpen(filename, NULL);
	if (!holdFire)
		launchBufferedReader(&readBigBed, data, &(data->bufferedReaderData));
}

void BigBedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigBedReaderData * data = (BigBedReaderData *) wi->data; 

	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}
	data->chrom = chrom;
	data->start = start;
	data->stop = finish;
	launchBufferedReader(&readBigBed, data, &(data->bufferedReaderData));
	wi->done = false;
	BigBedReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigBedReaderPop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
}

WiggleIterator * BigBedReader(char * f, bool holdFire) {
	BigBedReaderData * data = (BigBedReaderData *) calloc(1, sizeof(BigBedReaderData));
	openBigBed(data, f, holdFire);
	return newWiggleIterator(data, &BigBedReaderPop, &BigBedReaderSeek, 0, true);
}	
