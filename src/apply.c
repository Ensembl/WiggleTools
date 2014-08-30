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
#include <math.h>

#include "multiplexer.h"

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
	double default_value;
} BufferedWiggleIteratorData;

static BufferedWiggleIteratorData * createBufferedWiggleIteratorData(char * chrom, int start, int finish, float default_value) {
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
	bufferedData->default_value = default_value;
	return bufferedData;
}

void destroyBufferedWiggleIteratorData(BufferedWiggleIteratorData * data) {
	free(data->values);
	free(data->set);
	free(data);
}

void LooseBufferedWiggleIteratorPop(WiggleIterator * apply) {
	BufferedWiggleIteratorData * data = (BufferedWiggleIteratorData *) apply->data;
	if (apply->done)
		;
	else if (data->index == data->length)
		apply->done = true;
	else {
		apply->start = data->index;
		if (data->set[data->index]) {
			apply->value = data->values[data->index];
			data->index++;
		} else {
			apply->value = apply->default_value;
			while (data->index < data->length && !data->set[data->index])
				data->index++;
		}
		apply->finish = data->index;
	}
}

void StrictBufferedWiggleIteratorPop(WiggleIterator * apply) {
	BufferedWiggleIteratorData * data = (BufferedWiggleIteratorData *) apply->data;
	while (data->index < data->length) {
		if (data->set[data->index]) {
			apply->start = data->index;
			apply->finish = apply->start + 1;
			apply->value = data->values[data->index];
			data->index++;
			return;
		} else
			data->index++;
	}
	apply->done = true;
}

void BufferedWiggleIteratorSeek(WiggleIterator * apply, const char * chrom, int start, int finish) {
	fprintf(stderr, "Cannot seek on buffered iterator!");
	exit(1);
}

WiggleIterator * BufferedWiggleIterator(BufferedWiggleIteratorData * data, bool strict) {
	WiggleIterator * apply;
	if (strict)
		apply = newWiggleIterator(data, &StrictBufferedWiggleIteratorPop, &BufferedWiggleIteratorSeek, data->default_value);
	else
		apply = newWiggleIterator(data, &LooseBufferedWiggleIteratorPop, &BufferedWiggleIteratorSeek, data->default_value);
	apply->chrom = data->chrom;
	return apply;
}

//////////////////////////////////////////////////////
// Apply operator
//////////////////////////////////////////////////////

typedef struct applyWiggleIteratorData_st {
	WiggleIterator * regions;
	WiggleIterator * (**statistics)(WiggleIterator *);
	int profile_width;
	bool strict;
	WiggleIterator * input;
	BufferedWiggleIteratorData * head;
	BufferedWiggleIteratorData * tail;
} ApplyMultiplexerData;

static void createTarget(ApplyMultiplexerData * data) {
	BufferedWiggleIteratorData * bufferedData = createBufferedWiggleIteratorData(data->regions->chrom, data->regions->start, data->regions->finish, data->input->default_value);
	if (!data->head)
		data->head = bufferedData;
	else
		data->tail->next = bufferedData;
	data->tail = bufferedData;
}

static void createTargets(ApplyMultiplexerData * data) {
	// NOTE: the 10kb added allows the system to pull neighbouring regions together 
	while(!data->regions->done && (!data->head || (data->regions->start <= data->tail->finish + 1000000 && !strcmp(data->regions->chrom, data->tail->chrom)))) {
		if (data->regions->value)
			createTarget(data);
		pop(data->regions);
	}
}

static void pushDataOnBuffer(ApplyMultiplexerData * data, BufferedWiggleIteratorData * bufferedData) {
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

static void pushData(ApplyMultiplexerData * data) {
	BufferedWiggleIteratorData * bufferedData;

	for (bufferedData = data->head; bufferedData; bufferedData = bufferedData->next) {
		if (bufferedData->start >= data->input->finish)
			break;
		else
			pushDataOnBuffer(data, bufferedData);
	}
}

void ApplyMultiplexerPop(Multiplexer * apply) {
	ApplyMultiplexerData * data = (ApplyMultiplexerData *) apply->data;

	// If no ongoing jobs, create some
	if (data->head == NULL) {
		createTargets(data);
		if (data->head) {
			seek(data->input, data->head->chrom, data->head->start, data->tail->finish);	
		} else {
			apply->done = true;
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
	apply->chrom = data->head->chrom;
	apply->start = data->head->start;
	apply->finish = data->head->finish;
	WiggleIterator * wi = BufferedWiggleIterator(data->head, data->strict);
	if (data->statistics) {
		int i;
		for (i = apply->count-1; i >= 0; i--)
			wi = (data->statistics[i])(wi);
		runWiggleIterator(wi);
		i=0;
		while (wi->append) {
			apply->values[i] = *((double*) (wi->data));
			WiggleIterator * tmp = wi;
			wi = wi->append;
			if (tmp->data)
				free(tmp->data);
			free(tmp);
			i++;
		}
	} else
		regionProfile(wi, apply->values, apply->count, apply->finish - apply->start, false);

	// Discard struct
	BufferedWiggleIteratorData * bufferedData = data->head;
	if (data->tail == data->head) 
		data->tail = data->head = NULL;
	else
		data->head = data->head->next;
	destroyBufferedWiggleIteratorData(bufferedData);
	free(wi);
}

void ApplyMultiplexerSeek(Multiplexer * apply, const char * chrom, int start, int finish) {
	ApplyMultiplexerData * data = (ApplyMultiplexerData *) apply->data;
	BufferedWiggleIteratorData * bufferedData;
	while (data->head) {
		bufferedData = data->head;
		data->head = data->head->next;
		destroyBufferedWiggleIteratorData(bufferedData);
	}
	data->tail = NULL;
	seek(data->regions, chrom, start, finish);
}

Multiplexer * ApplyMultiplexer(WiggleIterator * regions, WiggleIterator * (**statistics)(WiggleIterator *), int count, WiggleIterator * dataset, bool strict) {
	ApplyMultiplexerData * data = (ApplyMultiplexerData *) calloc(1, sizeof(ApplyMultiplexerData));
	data->regions = regions;
	data->statistics = statistics;
	data->input = dataset;
	data->strict = strict;
	Multiplexer * res = newCoreMultiplexer(data, count, &ApplyMultiplexerPop, &ApplyMultiplexerSeek);
	int i;
	for (i=0; i < count; i++)
		res->default_values[i] = NAN;
	popMultiplexer(res);
	return res;
}

Multiplexer * ProfileMultiplexer(WiggleIterator * regions, int width, WiggleIterator * dataset) {
	ApplyMultiplexerData * data = (ApplyMultiplexerData *) calloc(1, sizeof(ApplyMultiplexerData));
	data->regions = regions;
	data->input = dataset;
	data->strict = false;
	Multiplexer * res = newCoreMultiplexer(data, width, &ApplyMultiplexerPop, &ApplyMultiplexerSeek);
	int i;
	for (i=0; i < width; i++)
		res->default_values[i] = NAN;
	popMultiplexer(res);
	return res;
}
