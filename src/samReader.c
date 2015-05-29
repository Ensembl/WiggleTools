// Copyright 2013-2014 EMBL-EBI
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

#include <stdlib.h>
#include <string.h> 

#include "wiggleIterator.h"
#include "fib.h"

typedef struct samReaderData_st {
	char  *filename;
	FILE * file;
	char * chrom;
	int stop;

	FibHeap * starts, * ends;
	char chromBuf[1000];
	char cigar[1000];
	int pos;
	bool done;
} SamReaderData;

static int isSeparator(char c) {
	char * separators = "MIDNSHPX=";
	int i;

	for (i = 0; i < 9; i++)
		if (c == separators[i])
			return 1;
	return 0;
}

static char * readNextCigarBlock(char * cigar, char * block) {
	int i;
	int lgth = strlen(cigar);
	for (i = 0; i < lgth && (i == 0 || !isSeparator(cigar[i-1])); i++) 
		block[i] = cigar[i];
	block[i] = '\0';

	if (i)
		return cigar + i;
	else
		return NULL;
}

static int storeReadComponent(FibHeap * starts, FibHeap * ends, int start, char * block) {
	int last = strlen(block) - 1;
	char type = block[last];
	block[last] = '\0';
	int count = atoi(block);
	
	switch (type) {
		case 'M':
		case 'X':
		case '=':
		case 'D':
			fh_insert(starts, start);
			fh_insert(ends, start + count);
		case 'N':
			return start + count;
		default:
			return start;
	}
}

static void readLine(WiggleIterator * wi) {
	SamReaderData * data = (SamReaderData *) wi->data;
	char line[5000];
	char chrom[1000];
	int pos;

	while (fgets(line, 5000, data->file)) {
		if (line[0] != '#' && line[0] != EOF) {
			sscanf(line, "%*s\t%*i\t%s\t%i\t%*i\t%s", chrom, &pos, data->cigar);

			if (strcmp(chrom, wi->chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && pos < data->pos)) {
				fprintf(stderr, "Sam file %s is not sorted!\nPosition %s:%i is before %s:%i\n", data->filename, chrom, pos, wi->chrom, data->pos);
				exit(1);
			}

			strcpy(data->chromBuf, chrom);
			data->pos = pos;
			return;
		}
	}

	data->done = true;
}

static void storeReadComponents(FibHeap * starts, FibHeap * ends, int start, char * cigar) {
	char block[100];
	char * ptr;

	for (ptr = readNextCigarBlock(cigar, block); ptr; ptr = readNextCigarBlock(ptr, block)) 
		start = storeReadComponent(starts, ends, start, block);
}

static void loadNextReads(WiggleIterator * wi) {
	SamReaderData * data = (SamReaderData *) wi->data;

	if (data->chromBuf[0] == '\0')
		readLine(wi);

	if (fh_empty(data->ends) && strcmp(wi->chrom, data->chromBuf)) {
		wi->chrom = (char *) calloc(strlen(data->chromBuf) + 1, sizeof(char));
		strcpy(wi->chrom, data->chromBuf);
	}

	while (!data->done && !strcmp(wi->chrom, data->chromBuf) && (fh_empty(data->ends) || fh_empty(data->starts) || data->pos <= fh_min(data->ends) || data->pos <= fh_min(data->starts))) {
		storeReadComponents(data->starts, data->ends, data->pos, data->cigar);
		readLine(wi);
	}
}

void SamReaderPop(WiggleIterator * wi) {
	SamReaderData * data = (SamReaderData *) wi->data;

	if (wi->done)
		return;

	if (wi->chrom[0] == '\0')
		wi->value = 0;

	while (fh_notempty(data->ends) && fh_min(data->ends) == wi->finish) {
		wi->value--;
		fh_extractmin(data->ends);
	}

	if (wi->value < 0) {
		fprintf(stderr, "Negative coverage at %s:%i???\n", wi->chrom, wi->finish);
		exit(1);
	}

	loadNextReads(wi);

	if (fh_notempty(data->ends)) {
		if (wi->value)
			wi->start = wi->finish;
		else 
			wi->start = fh_min(data->starts);

		while (fh_notempty(data->starts) && fh_min(data->starts) == wi->start) {
			wi->value++;
			fh_extractmin(data->starts);
		}

		if (fh_notempty(data->starts)) 
			wi->finish = fh_min(data->starts);
		if (fh_notempty(data->ends) && (fh_empty(data->starts) || fh_min(data->ends) < wi->finish)) 
			wi->finish = fh_min(data->ends);

		if (data->stop > 0) {
			if ((wi->start >= data->stop && strcmp(wi->chrom, data->chrom) == 0) || strcmp(wi->chrom, data->chrom) > 0)
				wi->done = true;
			else if (wi->finish > data->stop)
				wi->finish = data->stop;
		}

		if (wi->start < 0)
			abort();
		if (wi->start == 0)
			abort();
	}
	else {
		fclose(data->file);
		data->file = NULL;
		wi->done = true;
	}
}

void SamReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	SamReaderData * data = (SamReaderData*) wi->data;

	if (data->file == stdin) {
		fprintf(stderr, "Cannot do a seek on stdin stream!\n");
		exit(1);
	}

	data->stop = finish;
	data->chrom = chrom;

	if (!data->file || strcmp(chrom, wi->chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && start < wi->start)) {
		if (data->file)
			fclose(data->file);
		if (!(data->file = fopen(data->filename, "r"))) {
			fprintf(stderr, "Could not open input file %s\n", data->filename);
			exit(1);
		}
		wi->done = false;
		pop(wi);
	}

	fh_deleteheap(data->starts);
	fh_deleteheap(data->ends);
	data->starts = fh_makeheap();
	data->ends = fh_makeheap();
	data->done = true;
	wi->value = 0;
	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish < start))) 
		pop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
}

WiggleIterator * SamReader(char * filename) {
	SamReaderData * data = (SamReaderData *) calloc(1, sizeof(SamReaderData));
	data->filename = filename;
	data->stop = -1;
	data->starts = fh_makeheap();
	data->ends = fh_makeheap();
	if (strcmp(filename, "-")) {
		if (!(data->file = fopen(filename, "r"))) {
			fprintf(stderr, "Could not open input file %s\n", filename);
			exit(1);
		}
	} else
		data->file = stdin;
	return newWiggleIterator(data, &SamReaderPop, &SamReaderSeek, 0);
}
