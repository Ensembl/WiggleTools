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
	bool * set;
	struct bufferedWiggleIteratorData_st * next;
} BufferedWiggleIteratorData;

static BufferedWiggleIteratorData * createBufferedWiggleIteratorData(char * chrom, int start, int finish) {
	BufferedWiggleIteratorData * bufferedData = (BufferedWiggleIteratorData *) calloc(1, sizeof(BufferedWiggleIteratorData));
	if (!bufferedData) {
		fprintf(stderr, "Could not calloc %li bytes\n", sizeof(BufferedWiggleIteratorData));
		abort();
	}
	bufferedData->chrom = chrom;
	bufferedData->start = start;
	bufferedData->finish = finish;
	bufferedData->index = 0;
	bufferedData->length = finish - start;
	bufferedData->values = (double *) calloc(bufferedData->length, sizeof(double));
	bufferedData->set = (bool *) calloc(bufferedData->length, sizeof(bool));
	return bufferedData;
}

void destroyBufferedWiggleIteratorData(BufferedWiggleIteratorData * data) {
	free(data->values);
	free(data->set);
	free(data);
}

void BufferedWiggleIteratorPop(WiggleIterator * wi) {
	BufferedWiggleIteratorData * data = (BufferedWiggleIteratorData *) wi->data;
	while (++data->index < data->length) {
		if (data->set[data->index]) {
			wi->start = data->index;
			wi->finish = wi->start + 1;
			wi->value = data->values[data->index];
			return;
		}
	}
	wi->done = true;
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
	double * valuePtr;
	int profile_width;
	WiggleIterator * input;
	BufferedWiggleIteratorData * head;
	BufferedWiggleIteratorData * tail;
} ApplyWiggleIteratorData;

static void createTarget(ApplyWiggleIteratorData * data) {
	BufferedWiggleIteratorData * bufferedData = createBufferedWiggleIteratorData(data->regions->chrom, data->regions->start, data->regions->finish);
	if (!data->head)
		data->head = bufferedData;
	else
		data->tail->next = bufferedData;
	data->tail = bufferedData;
}

static void createTargets(ApplyWiggleIteratorData * data) {
	// NOTE: the 10kb added allows the system to pull neighbouring regions together 
	while(!data->regions->done && (!data->head || (data->regions->start <= data->tail->finish + 1000000 && !strcmp(data->regions->chrom, data->tail->chrom)))) {
		if (data->regions->value)
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

	for (pos = start; pos < finish; pos++) {
		bufferedData->values[pos] = data->input->value;
		bufferedData->set[pos] = true;
	}

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
	if (data->statistic)
		wi->value = data->statistic(BufferedWiggleIterator(data->head));
	else {
		wi->valuePtr = data->valuePtr;
		regionProfile(BufferedWiggleIterator(data->head), wi->valuePtr, data->profile_width, wi->finish - wi->start, false);
	}

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

WiggleIterator * ProfileWiggleIterator(WiggleIterator * regions, int width, WiggleIterator * dataset) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) calloc(1, sizeof(ApplyWiggleIteratorData));
	data->regions = regions;
	data->profile_width = width;
	data->input = dataset;
	data->valuePtr = calloc(width, sizeof(double));
	return newWiggleIterator(data, &ApplyWiggleIteratorPop, NULL);
}
