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

#include "wiggleMultiplexer.h"

typedef struct wiggleReducerData_st {
	Multiplexer * multi;
} WiggleReducerData;

void WiggleReducerSeek(WiggleIterator * iter, const char * chrom, int start, int finish) {
	WiggleReducerData * data = (WiggleReducerData* ) iter->data;
	seekMultiplexer(data->multi, chrom, start, finish);
}

////////////////////////////////////////////////////////
// Max
////////////////////////////////////////////////////////

void MaxReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	if (multi->inplay[0])
		wi->nextValue = multi->iters[0]->value;
	else
		wi->nextValue = 0;

	for (i = 1; i < multi->count; i++) {
		if (multi->inplay[i]) {
			if (wi->nextValue < multi->values[i])
				wi->nextValue = multi->values[i];
		} else {
			if (wi->nextValue < 0)
				wi->nextValue = 0;
		}
	}
	popMultiplexer(multi);
}

WiggleIterator * MaxReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MaxReductionPop, &WiggleReducerSeek);
}

WiggleIterator * MaxWiggleReducer(WiggleIterator** iters, int count) {
	return MaxReduction(newMultiplexer(iters, count));
}

////////////////////////////////////////////////////////
// Min
////////////////////////////////////////////////////////

void MinReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	if (multi->inplay[0])
		wi->nextValue = multi->iters[0]->value;
	else
		wi->nextValue = 0;

	for (i = 1; i < multi->count; i++) {
		if (multi->inplay[i]) {
			if (wi->nextValue > multi->values[i])
				wi->nextValue = multi->values[i];
		} else {
			if (wi->nextValue > 0)
				wi->nextValue = 0;
		}
	}
	popMultiplexer(multi);
}

static WiggleIterator * MinReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MinReductionPop, &WiggleReducerSeek);
}

WiggleIterator * MinWiggleReducer(WiggleIterator** iters, int count) {
	return MinReduction(newMultiplexer(iters, count));
}

////////////////////////////////////////////////////////
// Sum
////////////////////////////////////////////////////////

void SumReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	wi->nextValue = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			wi->nextValue += multi->values[i];
	popMultiplexer(multi);
}

static WiggleIterator * SumReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &SumReductionPop, &WiggleReducerSeek);
}

WiggleIterator * SumWiggleReducer(WiggleIterator** iters, int count) {
	return SumReduction(newMultiplexer(iters, count));
}

////////////////////////////////////////////////////////
// Product
////////////////////////////////////////////////////////

void ProductReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	wi->nextValue = 1;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i])
			wi->nextValue *= multi->values[i];
		else {
			wi->nextValue = 0;
			break;
		}
	}
	popMultiplexer(multi);
}

static WiggleIterator * ProductReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &ProductReductionPop, &WiggleReducerSeek);
}

WiggleIterator * ProductWiggleReducer(WiggleIterator** iters, int count) {
	return ProductReduction(newMultiplexer(iters, count));
}

////////////////////////////////////////////////////////
// Mean
////////////////////////////////////////////////////////

void MeanReductionPop(WiggleIterator * wi) {
	int i;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	wi->nextValue = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			wi->nextValue += multi->values[i];
	wi->nextValue /= multi->count;
	popMultiplexer(multi);
}

static WiggleIterator * MeanReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &MeanReductionPop, &WiggleReducerSeek);
}

WiggleIterator * MeanWiggleReducer(WiggleIterator** iters, int count) {
	return MeanReduction(newMultiplexer(iters, count));
}


////////////////////////////////////////////////////////
// Variance
////////////////////////////////////////////////////////

void VarianceReductionPop(WiggleIterator * wi) {
	int i;
	double mean, diff;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	mean = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			mean += multi->values[i];
	mean /= multi->count;
	
	wi->nextValue = 0;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i]) {
			diff = (mean - multi->values[i]);
			wi->nextValue += diff * diff;
		} else {
			wi->nextValue += mean * mean;
		}
	}
	wi->nextValue /= multi->count;
	
	popMultiplexer(multi);
}

static WiggleIterator * VarianceReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &VarianceReductionPop, &WiggleReducerSeek);
}

WiggleIterator * VarianceWiggleReducer(WiggleIterator** iters, int count) {
	return VarianceReduction(newMultiplexer(iters, count));
}

////////////////////////////////////////////////////////
// StdDev
////////////////////////////////////////////////////////

void StdDevReductionPop(WiggleIterator * wi) {
	int i;
	double mean, diff;

	if (wi->nextDone)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	mean = 0;
	for (i = 0; i < multi->count; i++)
		if (multi->inplay[i])
			mean += multi->values[i];
	mean /= multi->count;
	
	wi->nextValue = 0;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i]) {
			diff = (mean - multi->values[i]);
			wi->nextValue += diff * diff;
		} else {
			wi->nextValue += mean * mean;
		}
	}
	wi->nextValue /= multi->count;
	wi->nextValue = sqrt(wi->nextValue);
	
	popMultiplexer(multi);
}

static WiggleIterator * StdDevReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	return newWiggleIterator(data, &StdDevReductionPop, &WiggleReducerSeek);
}

WiggleIterator * StdDevWiggleReducer(WiggleIterator** iters, int count) {
	return StdDevReduction(newMultiplexer(iters, count));
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

	if (wi->nextDone)
		return;

	MedianWiggleReducerData * data = (MedianWiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->nextDone = true;
		return;
	}

	wi->nextChrom = multi->chrom;
	wi->nextStart = multi->start;
	wi->nextFinish = multi->finish;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i])
			data->vals[i] = multi->values[i];
		else
			data->vals[i] = 0;
	}
	qsort(data->vals, multi->count, sizeof(double), &compDoubles);
	wi->nextValue = data->vals[multi->count / 2];
	
	popMultiplexer(multi);
}

static WiggleIterator * MedianReduction(Multiplexer * multi) {
	MedianWiggleReducerData * data = (MedianWiggleReducerData *) calloc(1, sizeof(MedianWiggleReducerData));
	data->multi = multi;
	data->vals = (double *) calloc(data->multi->count, sizeof(double));
	return newWiggleIterator(data, &MedianReductionPop, &MedianWiggleReducerSeek);
}

WiggleIterator * MedianWiggleReducer(WiggleIterator** iters, int count) {
	return MedianReduction(newMultiplexer(iters, count));
}