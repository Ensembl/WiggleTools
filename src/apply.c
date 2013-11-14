// Copyright (c) 2013, Daniel Zerbino
// All rights valueerved.
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

// Local header
#include <stdlib.h>
#include <string.h>

#include "wiggleIterator.h"

//////////////////////////////////////////////////////
// Buffered wiggleIterator
//////////////////////////////////////////////////////

typedef struct bufferedWiggleIteratorData_st {
	char * chrom;
	int start;
	int finish;
	int index;
	int length;
	double * values;
	struct bufferedWiggleIteratorData_st * next;
} BufferedWiggleIteratorData;

void destroyBufferedWiggleIteratorData(BufferedWiggleIteratorData * data) {
	free(data->values);
	free(data);
}

void BufferedWiggleIteratorPop(WiggleIterator * wi) {
	BufferedWiggleIteratorData * data = (BufferedWiggleIteratorData *) wi->data;
	if (data->index++ < data->length) {
		wi->start = data->start + data->index;
		wi->finish = wi->start + 1;
		wi->value = data->values[data->index];
	} else {
		wi->done = true;
	}
}

void BufferedWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	fprintf(stderr, "Cannot seek on buffered iterator!");
	exit(1);
}

WiggleIterator * BufferedWiggleIterator(BufferedWiggleIteratorData * data) {
	WiggleIterator * wi = newWiggleIterator(data, &BufferedWiggleIteratorPop, &BufferedWiggleIteratorSeek);
	wi->chrom = data->chrom;
	return wi;
}

//////////////////////////////////////////////////////
// Apply operator
//////////////////////////////////////////////////////

typedef struct applyWiggleIteratorData_st {
	WiggleIterator * regions;
	double (*statistic)(WiggleIterator *);
	WiggleIterator * input;
	BufferedWiggleIteratorData * head;
	BufferedWiggleIteratorData * tail;
} ApplyWiggleIteratorData;

static void createTarget(ApplyWiggleIteratorData * data) {
	BufferedWiggleIteratorData * bufferedData = (BufferedWiggleIteratorData *) calloc(1, sizeof(BufferedWiggleIteratorData));
	if (!bufferedData) {
		fprintf(stderr, "Could not calloc %li bytes\n", sizeof(BufferedWiggleIteratorData));
		abort();
	}
	bufferedData->chrom = data->regions->chrom;
	bufferedData->start = data->regions->start;
	bufferedData->finish = data->regions->finish;
	bufferedData->index = 0;
	bufferedData->length = data->regions->finish - data->regions->start;
	bufferedData->values = (double *) calloc(bufferedData->length, sizeof(double));
	if (!data->head)
		data->head = bufferedData;
	else
		data->tail->next = bufferedData;
	data->tail = bufferedData;
}

static void createTargets(ApplyWiggleIteratorData * data) {
	while(!data->regions->done && (!data->head || (data->regions->start <= data->tail->finish && !strcmp(data->regions->chrom, data->tail->chrom)))) {
		createTarget(data);
		pop(data->regions);
	}
}

static void pushDataOnBuffer(ApplyWiggleIteratorData * data, BufferedWiggleIteratorData * bufferedData) {
	int pos, start, finish;

	if (bufferedData->start > data->input->start)
		start = 0;
	else
		start = data->input->start - bufferedData->start;

	if (bufferedData->finish < data->input->finish)
		finish = bufferedData->length;
	else
		finish = data->input->finish - bufferedData->start;

	for (pos = start; pos < finish; pos++)
		bufferedData->values[pos] = data->input->value;

}

static void pushData(ApplyWiggleIteratorData * data) {
	BufferedWiggleIteratorData * bufferedData;

	for (bufferedData = data->head; bufferedData; bufferedData = bufferedData->next) {
		if (bufferedData->start >= data->input->finish)
			break;
		else
			pushDataOnBuffer(data, bufferedData);
	
	}

}

void ApplyWiggleIteratorPop(WiggleIterator * wi) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) wi->data;

	// If no ongoing jobs, create some
	if (data->head == NULL) {
		createTargets(data);
		if (data->head) {
			seek(data->input, data->head->chrom, data->head->start, data->tail->finish);	
		} else {
			wi->done = true;
			return;
		}
	}

	// If ongoing jobs are running:
	// Push enough data to finish the first job
	while (!data->input->done && !strcmp(data->input->chrom, data->head->chrom) && data->input->start < data->head->finish) {
		pushData(data);
		pop(data->input);
	}

	// Return value
	wi->chrom = data->head->chrom;
	wi->start = data->head->start;
	wi->finish = data->head->finish;
	wi->value = data->statistic(BufferedWiggleIterator(data->head));

	// Discard struct
	BufferedWiggleIteratorData * bufferedData = data->head;
	if (data->tail == data->head) 
		data->tail = data->head = NULL;
	else
		data->head = data->head->next;
	destroyBufferedWiggleIteratorData(bufferedData);
}

void ApplyWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) wi->data;
	BufferedWiggleIteratorData * bufferedData;
	while (data->head) {
		bufferedData = data->head;
		data->head = data->head->next;
		destroyBufferedWiggleIteratorData(bufferedData);
	}
	data->tail = NULL;
	seek(data->regions, chrom, start, finish);
}

WiggleIterator * ApplyWiggleIterator(WiggleIterator * regions, double (*statistic)(WiggleIterator *), WiggleIterator * dataset) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) calloc(1, sizeof(ApplyWiggleIteratorData));
	data->regions = regions;
	data->statistic = statistic;
	data->input = dataset;
	return newWiggleIterator(data, &ApplyWiggleIteratorPop, NULL);
}
