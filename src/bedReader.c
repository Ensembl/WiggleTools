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

#include <stdlib.h>
#include <string.h> 

#include "wiggleIterators.h"

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

	if (wi->done)
		return;

	if (fgets(line, 5000, data->file)) {
		sscanf(line, "%s\t%i\t%i", chrom, &wi->start, &wi->finish);
		// I like my finishes to be non inclusive...
		wi->finish++;

		// The reason for creating a new string instead of simply 
		// overwriting is that other functions may still be pointin
		// at the old label
		if (wi->chrom[0] == '\0' || strcmp(wi->chrom, chrom)) {
			wi->chrom = (char *) calloc(strlen(chrom), sizeof(char));
			strcpy(wi->chrom, chrom);
		}

		if (data->stop > 0 && (strcmp(wi->chrom, data->chrom) > 0 || (strcmp(data->chrom, wi->chrom) == 0 && wi->start > data->stop))) {
			wi->done = true;
		}
	} else {
		fclose(data->file);
		data->file = NULL;
		wi->done = true;
	}
}

void BedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BedReaderData * data = (BedReaderData*) wi->data;

	if (wi->done || strcmp(chrom, wi->chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && start < wi->start)) {
		if (data->file)
			fclose(data->file);
		if (!(data->file = fopen(data->filename, "r"))) {
			printf("Could not open input file %s\n", data->filename);
			exit(1);
		}
		pop(wi);
	}

	data->chrom = chrom;
	data->stop = finish;
	wi->done = false;
	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish < start))) 
		pop(wi);
}

WiggleIterator * BedReader(char * filename) {
	BedReaderData * data = (BedReaderData *) calloc(1, sizeof(BedReaderData));
	data->filename = filename;
	data->stop = -1;
	if (!(data->file = fopen(filename, "r"))) {
		printf("Could not open bed file %s\n", filename);
		exit(1);
	}
	return UnionWiggleIterator(newWiggleIterator(data, &BedReaderPop, &BedReaderSeek));
}
