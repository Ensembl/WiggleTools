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

#include <stdlib.h>
#include <string.h> 

#include "wiggleIterator.h"

typedef struct bedReaderData_st {
	char  *filename;
	FILE * file;
	char * chrom;
	int stop;
} BedReaderData;

void BedReaderPop(WiggleIterator * wi) {
	BedReaderData * data = (BedReaderData *) wi->data;
	char line[5000];
	char chrom[1000];
	char sign;

	if (wi->done)
		return;

	if (fgets(line, 5000, data->file)) {
		sscanf(line, "%s\t%i\t%i\t%*s\t%c", chrom, &wi->start, &wi->finish, &sign);

		if (sign == '+')
			wi->strand = 1;
		else if (sign == '-')
			wi->strand = -1;
		else
			wi->strand = 0;

		// Conversion from 0 to 1-based...
		wi->start++;
		wi->finish++;

		// The reason for creating a new string instead of simply 
		// overwriting is that other functions may still be pointin
		// at the old label
		if (wi->chrom[0] == '\0' || strcmp(wi->chrom, chrom)) {
			wi->chrom = (char *) calloc(strlen(chrom), sizeof(char));
			strcpy(wi->chrom, chrom);
		}

		if (data->stop > 0) {
			if ((wi->start >= data->stop && strcmp(wi->chrom, data->chrom) == 0) || strcmp(wi->chrom, data->chrom) > 0) {
				wi->done = true;
				return;
			} else if (wi->finish > data->stop) {
				wi->finish = data->stop;
			}
		}
	} else {
		fclose(data->file);
		data->file = NULL;
		wi->done = true;
	}
}

void BedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BedReaderData * data = (BedReaderData*) wi->data;

	data->stop = finish;
	data->chrom = chrom;

	if (wi->done || strcmp(chrom, wi->chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && start < wi->start)) {
		if (data->file)
			fclose(data->file);
		if (!(data->file = fopen(data->filename, "r"))) {
			fprintf(stderr, "Could not open input file %s\n", data->filename);
			exit(1);
		}
		wi->done = false;
		pop(wi);
	}

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish < start))) 
		pop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
}

WiggleIterator * BedReader(char * filename) {
	BedReaderData * data = (BedReaderData *) calloc(1, sizeof(BedReaderData));
	data->filename = filename;
	data->stop = -1;
	if (!(data->file = fopen(filename, "r"))) {
		fprintf(stderr, "Could not open bed file %s\n", filename);
		exit(1);
	}
	WiggleIterator * res = newWiggleIterator(data, &BedReaderPop, &BedReaderSeek);
	res->overlaps = true;
	return res;
}
