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

#include "wiggleTools.h"
#include "wiggleIterators.h"

//////////////////////////////////////////////////////
// File Reader
//////////////////////////////////////////////////////

typedef struct wiggleReaderData_st {
	FILE * infile;
	size_t length;
	char * line;
	bool fixedStep;
	int step;
	int span;
} WiggleReaderData;

static void WiggleReaderReadHeader(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData *) wi->data;
	bool chrom_b = true;
	bool start_b = true;
	bool step_b = true;
	const char * seps = " \t=";
	char * token = strtok(data->line, seps);
	if (!strcmp(token, "variableStep"))
		data->fixedStep = false;
	else if (!strcmp(token, "fixedStep"))
		data->fixedStep = true;
	else {
		printf("Header line wi->starting with word %s\n", token);
		exit(1);
	}

	// Default
	data->span = 1;

	// Reading following parameters
	token = strtok(NULL, seps);
	while(token) {
		if (!strcmp(token, "chrom")) {
			chrom_b = false;
			token = strtok(NULL, seps);
			if (!token) {
				printf("Empty wi->chromosome name!\n");
				exit(1);
			}
			strcpy(wi->chrom, token);
		}
		if (!strcmp(token, "start")) {
			start_b = false;
			token = strtok(NULL, seps);
			if (!token) {
				printf("Empty wi->start position!\n");
				exit(1);
			}
			sscanf(token, "%i", &(wi->start));
		}
		if (!strcmp(token, "span")) {
			token = strtok(NULL, seps);
			if (!token) {
				printf("Empty span length!\n");
				exit(1);
			}
			sscanf(token, "%i", &(data->span));
		}
		if (!strcmp(token, "step")) {
			step_b = false;
			if (!data->fixedStep) {
				printf("Cannot specify step length on a variable length track\n");
				exit(1);
			}
			token = strtok(NULL, seps);
			if (!token) {
				printf("Empty step length!\n");
				exit(1);
			}
			sscanf(token, "%i", &(data->step));
		}
		token = strtok(NULL, seps);
	}

	// Checking that all compulsory fields were filled:
	if ((data->fixedStep && (chrom_b || start_b || step_b)) || (!data->fixedStep && chrom_b)) {
		printf("Invalid header, missing data: %s\n", data->line);
		exit(1);
	}

	// Backing off so as not to offset the first line
	if (data->fixedStep)
		wi->start -= data->step;
}

static void WiggleReaderReadFixedStepLine(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	sscanf(data->line, "%lf", &(wi->value));
	wi->start += data->step;
	wi->finish = wi->start + data->span;
}

static void WiggleReaderReadVariableStepLine(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	sscanf(data->line, "%i\t%lf", &(wi->start), &(wi->value));
	wi->finish = wi->start + data->span;
}

static void WiggleReaderReadLine(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	if (data->fixedStep)
		WiggleReaderReadFixedStepLine(wi);
	else
		WiggleReaderReadVariableStepLine(wi);
}

static bool WiggleReaderPopLine(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	if (fgets(data->line, data->length, data->infile)) {
		return true;
	} else {
		wi->done = true;
		fclose(data->infile);
		return false;
	}
}

static void WiggleReaderPop(WiggleIterator * wi) {
	WiggleReaderData * data = (WiggleReaderData*) wi->data;
	while (WiggleReaderPopLine(wi)) { 
		if (data->line[0] == 'v' || data->line[0] == 'f') {
			WiggleReaderReadHeader(wi);
		} else {
			WiggleReaderReadLine(wi);
			break;
		}
	}
}

WiggleIterator * WiggleReader(char * f) {
	WiggleReaderData * data = (WiggleReaderData *) calloc(1, sizeof(WiggleReaderData));
	if (strcmp(f, "-"))
		data->infile = openOrFail(f, "input file", "r");
	else
		data->infile = stdin;
	data->length = 5000;
	data->line = (char *) calloc(data->length, 1);
	return newWiggleIterator(data, &WiggleReaderPop);
}	

