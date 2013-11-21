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
#include <stdio.h>
#include <string.h>

#include "wiggleIterator.h"

//////////////////////////////////////////////////////
// File Reader
//////////////////////////////////////////////////////

enum readingMode {FIXED_STEP, VARIABLE_STEP, BED_GRAPH};

typedef struct wiggleReaderData_st {
	char * filename;
	FILE * file;
	enum readingMode readingMode;
	int step;
	int span;
	char words[5];
	int stop;
} WiggleReaderData;


static void WiggleReaderReadHeader(WiggleIterator * wi, WiggleReaderData * data, char * line) {
	bool chrom_b = true;
	bool start_b = true;
	bool step_b = true;
	const char * seps = " \t=";
	char * token = strtok(line, seps);

	// Default
	data->span = 1;

	// Reading following parameters
	token = strtok(NULL, seps);
	while(token) {
		if (!strcmp(token, "chrom")) {
			chrom_b = false;
			token = strtok(NULL, seps);
			if (!token) {
				fprintf(stderr, "Empty wi->chromosome name!\n");
				exit(1);
			}
			strcpy(wi->chrom, token);
		}
		if (!strcmp(token, "start")) {
			start_b = false;
			token = strtok(NULL, seps);
			if (!token) {
				fprintf(stderr, "Empty wi->start position!\n");
				exit(1);
			}
			sscanf(token, "%i", &(wi->start));
		}
		if (!strcmp(token, "span")) {
			token = strtok(NULL, seps);
			if (!token) {
				fprintf(stderr, "Empty span length!\n");
				exit(1);
			}
			sscanf(token, "%i", &(data->span));
		}
		if (!strcmp(token, "step")) {
			step_b = false;
			if (data->readingMode == VARIABLE_STEP) {
				fprintf(stderr, "Cannot specify step length on a variable length track\n");
				exit(1);
			}
			token = strtok(NULL, seps);
			if (!token) {
				fprintf(stderr, "Empty step length!\n");
				exit(1);
			}
			sscanf(token, "%i", &(data->step));
		}
		token = strtok(NULL, seps);
	}

	// Checking that all compulsory fields were filled:
	if ((data->readingMode == FIXED_STEP && (chrom_b || start_b || step_b)) || (data->readingMode == VARIABLE_STEP && chrom_b)) {
		fprintf(stderr, "Invalid header, missing data: %s\n", line);
		exit(1);
	}

	// Backing off so as not to offset the first line
	if (data->readingMode == FIXED_STEP)
		wi->start -= data->step;
}

static void WiggleReaderReadFixedStepLine(WiggleIterator * wi, char * line, int step, int span) {
	sscanf(line, "%lf", &(wi->value));
	wi->start += step;
	wi->finish = wi->start + span;
}

static void WiggleReaderReadVariableStepLine(WiggleIterator * wi, char * line, int span) {
	sscanf(line, "%i\t%lf", &(wi->start), &(wi->value));
	wi->finish = wi->start + span;
}

static void WiggleReaderReadBedGraphLine(WiggleIterator * wi, char * line) {
	char * buffer = calloc(500, sizeof(char));
	sscanf(line, "%s\t%i\t%i\t%lf", buffer, &(wi->start), &(wi->finish), &(wi->value));
	// BedGraphs are 0 based, half open
	wi->start++;
	wi->finish++;
	if (strcmp(buffer, wi->chrom))
		wi->chrom = buffer;
	else
		free(buffer);
}

static int countWords(char * line) {
	int count = 0;
	char * ptr;

	if (line[0] != ' ' && line[0] != '\t')
		count++;

	for (ptr = line; *ptr; ptr++)
		if (*ptr == ' ' || *ptr == '\t')
			count++;

	return count;
}

static void WiggleReaderPop(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	char line[5000];

	if (wi->done)
		return;

	while (fgets(line, 5000, data->file)) {
		if ( !strncmp("variableStep", line, 12)) {
			data->readingMode = VARIABLE_STEP;
			WiggleReaderReadHeader(wi, data, line);
			continue;
		} else if (!strncmp("fixedStep", line, 9)) {
			data->readingMode = FIXED_STEP;
			WiggleReaderReadHeader(wi, data, line);
			continue;
		} 
		
		switch (countWords(line)) {
		case 4:
			data->readingMode = BED_GRAPH;
			WiggleReaderReadBedGraphLine(wi, line);
			break;
		case 2:
			if (data->readingMode != VARIABLE_STEP) {
				fprintf(stderr, "Badly formatted fixed step line:\n%s", line);
				exit(1);
			}
			WiggleReaderReadVariableStepLine(wi, line,data->span);
			break;
		case 1:
			if (data->readingMode != FIXED_STEP) {
				fprintf(stderr, "Badly formatted variable step line:\n%s", line);
				exit(1);
			}
			WiggleReaderReadFixedStepLine(wi, line, data->step, data->span);
			break;
		default:
			fprintf(stderr, "Badly formatted wiggle or bed graph line :\n%s", line);
			exit(1);

		}

		if (data->stop > 0 && wi->start > data->stop)
			wi->done = true;

		return;

	}
	fclose(data->file);
	data->file = NULL;
	wi->done = true;
}

void WiggleReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;

	wi->chrom = chrom;
	data->stop = finish;

	if (strcmp(chrom, wi->chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && start < wi->start)) {
		if (data->file)
			fclose(data->file);
		if (!(data->file = fopen(data->filename, "r"))) {
			printf("Cannot open input file %s\n", data->filename);
			exit(1);
		}
		WiggleReaderPop(wi);
	}

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		WiggleReaderPop(wi);

}

WiggleIterator * WiggleReader(char * f) {
	WiggleReaderData * data = (WiggleReaderData *) calloc(1, sizeof(WiggleReaderData));
	data->filename = f;
	if (strcmp(f, "-")) {
		if (!(data->file = fopen(f, "r"))) {
			printf("Could not open input file %s\n", f);
			exit(1);
		}
	} else
		data->file = stdin;
	data->readingMode = BED_GRAPH;
	data->stop = -1;
	return CompressionWiggleIterator(newWiggleIterator(data, &WiggleReaderPop, &WiggleReaderSeek));
}	
