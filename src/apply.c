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

// Local header
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "multiplexer.h"

const int MAX_BUFFER = 1e6;
const int MAX_BUFFER_SUM = 1e6;
const int MAX_SEEK = 10;

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
	if (data->values) {
		free(data->values);
		free(data->set);
	}
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
		apply = newWiggleIterator(data, &StrictBufferedWiggleIteratorPop, &BufferedWiggleIteratorSeek, data->default_value, false);
	else
		apply = newWiggleIterator(data, &LooseBufferedWiggleIteratorPop, &BufferedWiggleIteratorSeek, data->default_value, false);
	apply->chrom = data->chrom;
	return apply;
}

//////////////////////////////////////////////////////
// Fill in wiggleIterator
//////////////////////////////////////////////////////

typedef struct fillInUnaryData_st {
	char * chrom;
	int start;
	int finish;
	bool first;
	WiggleIterator * source;
} FillInUnaryData;

void FillInUnaryPop(WiggleIterator * wi) {
	FillInUnaryData * data = (FillInUnaryData*) wi->data;
	WiggleIterator * source = data->source;

	if (data->first) {
		wi->chrom = data->chrom;
		wi->finish = data->start;
		data->first = false;
	}
	
	if (source->done) {
		if (wi->finish == data->finish)
			wi->done = true;
		else {
			wi->start = wi->finish;
			wi->finish = data->finish;
			wi->value = wi->default_value;
		}
	} else if (source->start > wi->finish) {
		wi->start = wi->finish;
		wi->finish = source->start;
		wi->value = wi->default_value;
	} else {
		wi->start = source->start;
		wi->finish = source->finish;
		wi->value = source->value;
		pop(data->source);
	}
}

WiggleIterator * FillInUnaryWiggleIterator(WiggleIterator * source, char * chrom, int start, int finish) {
	FillInUnaryData * data = (FillInUnaryData*) calloc(1, sizeof(FillInUnaryData));
	data->chrom = chrom;
	data->start = start;
	data->finish = finish;
	data->source = source;
	data->first = true;
	seek(source, chrom, start, finish);
	return newWiggleIterator(data, &FillInUnaryPop, NULL, source->default_value, false);
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

static BufferedWiggleIteratorData * createTarget(ApplyMultiplexerData * data) {
	return createBufferedWiggleIteratorData(data->regions->chrom, data->regions->start, data->regions->finish, data->input->default_value);
}

static void addTarget(ApplyMultiplexerData * data, BufferedWiggleIteratorData * bufferedData) {
	if (!data->head)
		data->head = bufferedData;
	else
		data->tail->next = bufferedData;
	data->tail = bufferedData;
}

static void createTargets(ApplyMultiplexerData * data) {
	if (data->regions->finish - data->regions->start >= MAX_BUFFER) {
		addTarget(data, createTarget(data));
		pop(data->regions);
	} else {
		int length;
		int last_finish = data->regions->finish;
		int total_buffers = 0;
		while(!data->regions->done 
		      && (length = data->regions->finish - data->regions->start) < MAX_BUFFER
		      && (!data->head 
			  || ((total_buffers += length) < MAX_BUFFER_SUM && data->regions->start <= last_finish + MAX_SEEK && !strcmp(data->regions->chrom, data->tail->chrom))
			 )
		     ) 
		{
			if (data->regions->finish > last_finish)
				last_finish = data->regions->finish;
			addTarget(data, createTarget(data));
			pop(data->regions);
		}
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

BufferedWiggleIteratorData * popApplyMultiplexerData(ApplyMultiplexerData * data) {
	BufferedWiggleIteratorData * bufferedData = data->head;
	if (data->tail == data->head) 
		data->tail = data->head = NULL;
	else
		data->head = data->head->next;
	return bufferedData;
}

void computeApplyValues(Multiplexer * apply, ApplyMultiplexerData * data, BufferedWiggleIteratorData * bufferedData) {
	WiggleIterator * wi;
	if (bufferedData->values)
		wi = BufferedWiggleIterator(bufferedData, data->strict);
	else if (data->strict) {
		wi = data->input;
		seek(wi, bufferedData->chrom, bufferedData->start, bufferedData->finish);
	} else
		wi = FillInUnaryWiggleIterator(data->input, bufferedData->chrom, bufferedData->start, bufferedData->finish);

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
			destroyWiggleIterator(tmp);
			i++;
		}
	} else
		regionProfile(wi, apply->values, apply->count, apply->finish - apply->start, false);

	if (wi != data->input) {
		// Careful not to destroy buffered data. It requires special function and is destroyed elsewhere.
		if (wi->data != bufferedData)
			free(wi->data);
		free(wi);
	}
}

void  updateApplyMultiplexer(Multiplexer * apply, ApplyMultiplexerData * data, BufferedWiggleIteratorData * bufferedData) {
	apply->chrom = bufferedData->chrom;
	apply->start = bufferedData->start;
	apply->finish = bufferedData->finish;
	computeApplyValues(apply, data, bufferedData);
}

void ApplyMultiplexerPop(Multiplexer * apply) {
	ApplyMultiplexerData * data = (ApplyMultiplexerData *) apply->data;

	// If no ongoing jobs, create some
	if (data->head == NULL) {
		// Note: only exit if no more regions AND no targets waiting 
		if (data->regions->done) {
			apply->done = true;
			return;
		} 
		createTargets(data);
		seek(data->input, data->head->chrom, data->head->start, data->tail->finish);	
	}

	// If ongoing targets are reading:
	// Push enough data to finish the first job
	if (data->head->values) {
		while (!data->input->done && data->input->start < data->head->finish && !strcmp(data->input->chrom, data->head->chrom)) {
			pushData(data);
			pop(data->input);
		}
	}

	// Return value
	BufferedWiggleIteratorData * bufferedData = popApplyMultiplexerData(data);
	updateApplyMultiplexer(apply, data, bufferedData);
	destroyBufferedWiggleIteratorData(bufferedData);
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
