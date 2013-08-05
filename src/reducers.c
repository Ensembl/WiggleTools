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
#include <math.h>

#include "multiplexer.h"

typedef struct wiggleReducerData_st {
	Multiplexer * multi;
} WiggleReducerData;

void WiggleReducerSeek(WiggleIterator * iter, const char * chrom, int start, int finish) {
	WiggleReducerData * data = (WiggleReducerData* ) iter->data;
	seekMultiplexer(data->multi, chrom, start, finish);
	pop(iter);
}

////////////////////////////////////////////////////////
// Select
////////////////////////////////////////////////////////

typedef struct wiggleSelectData_st {
	Multiplexer * multi;
	int index;
} WiggleSelectData;

void SelectReductionPop(WiggleIterator * wi) {
	if (wi->done)
		return;

	WiggleSelectData * data = (WiggleSelectData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	while (!multi->inplay[data->index]) {
		popMultiplexer(multi);
		if (multi->done) {
			wi->done = true;
			return;
		}
	}

	wi->value = multi->iters[data->index]->value;
	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	popMultiplexer(multi);
}

WiggleIterator * SelectReduction(Multiplexer * multi, int index) {
	WiggleSelectData * data = (WiggleSelectData *) calloc(1, sizeof(WiggleSelectData));
	data->multi = multi;
	data->index = index;
	return newWiggleIterator(data, &SelectReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Max
////////////////////////////////////////////////////////

void MaxReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	if (multi->inplay[0])
		wi->value = multi->iters[0]->value;
	else
		wi->value = 0;

	for (i = 1; i < multi->count; i++) {
		if (multi->inplay[i]) {
			if (wi->value < multi->values[i])
				wi->value = multi->values[i];
		} else {
			if (wi->value < 0)
				wi->value = 0;
		}
	}
	popMultiplexer(multi);
}

WiggleIterator * MaxReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MaxReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Min
////////////////////////////////////////////////////////

void MinReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	if (multi->inplay[0])
		wi->value = multi->iters[0]->value;
	else
		wi->value = 0;

	for (i = 1; i < multi->count; i++) {
		if (multi->inplay[i]) {
			if (wi->value > multi->values[i])
				wi->value = multi->values[i];
		} else {
			if (wi->value > 0)
				wi->value = 0;
		}
	}
	popMultiplexer(multi);
}

WiggleIterator * MinReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MinReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Sum
////////////////////////////////////////////////////////

void SumReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	wi->value = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			wi->value += multi->values[i];
	popMultiplexer(multi);
}

WiggleIterator * SumReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &SumReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Product
////////////////////////////////////////////////////////

void ProductReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	wi->value = 1;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i])
			wi->value *= multi->values[i];
		else {
			wi->value = 0;
			break;
		}
	}
	popMultiplexer(multi);
}

WiggleIterator * ProductReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &ProductReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Mean
////////////////////////////////////////////////////////

void MeanReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	wi->value = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			wi->value += multi->values[i];
	wi->value /= multi->count;
	popMultiplexer(multi);
}

WiggleIterator * MeanReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MeanReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Variance
////////////////////////////////////////////////////////

void VarianceReductionPop(WiggleIterator * wi) {
	int i;
	double mean, diff;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	mean = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			mean += multi->values[i];
	mean /= multi->count;
	
	wi->value = 0;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i]) {
			diff = (mean - multi->values[i]);
			wi->value += diff * diff;
		} else {
			wi->value += mean * mean;
		}
	}
	wi->value /= multi->count;
	
	popMultiplexer(multi);
}

WiggleIterator * VarianceReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &VarianceReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// StdDev
////////////////////////////////////////////////////////

void StdDevReductionPop(WiggleIterator * wi) {
	int i;
	double mean, diff;

	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	mean = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			mean += multi->values[i];
	mean /= multi->count;
	
	wi->value = 0;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i]) {
			diff = (mean - multi->values[i]);
			wi->value += diff * diff;
		} else {
			wi->value += mean * mean;
		}
	}
	wi->value /= multi->count;
	wi->value = sqrt(wi->value);
	
	popMultiplexer(multi);
}

WiggleIterator * StdDevReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &StdDevReductionPop, &WiggleReducerSeek);
}

////////////////////////////////////////////////////////
// Median
////////////////////////////////////////////////////////

typedef struct medianMiggleReducerData_st {
	Multiplexer * multi;
	double * vals;
} MedianWiggleReducerData;

void MedianWiggleReducerSeek(WiggleIterator * iter, const char * chrom, int start, int finish) {
	MedianWiggleReducerData * data = (MedianWiggleReducerData* ) iter->data;
	seekMultiplexer(data->multi, chrom, start, finish);
}

static int compDoubles(const void * A, const void * B) {
	double * a = (double *) A;
	double * b = (double *) B;
	if (*a < *b)
		return -1;
	else if (*a > *b)
		return 1;
	else
		return 0;
}

void MedianReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->done)
		return;

	MedianWiggleReducerData * data = (MedianWiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i])
			data->vals[i] = multi->values[i];
		else
			data->vals[i] = 0;
	}
	qsort(data->vals, multi->count, sizeof(double), &compDoubles);
	wi->value = data->vals[multi->count / 2];
	
	popMultiplexer(multi);
}

WiggleIterator * MedianReduction(Multiplexer * multi) {
	MedianWiggleReducerData * data = (MedianWiggleReducerData *) calloc(1, sizeof(MedianWiggleReducerData));
	data->multi = multi;
	data->vals = (double *) calloc(data->multi->count, sizeof(double));
	return newWiggleIterator(data, &MedianReductionPop, &MedianWiggleReducerSeek);
}
