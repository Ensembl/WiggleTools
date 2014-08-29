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

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Local header
#include "wiggleIterator.h"
#include "multiplexer.h"

//////////////////////////////////////////////////////
// Generic function for all statistics
//
// Note: If creating alternate statistic with 
// richer data structs, remember to always start
// with the result stored as double.
//////////////////////////////////////////////////////

typedef struct statData {
	double res;
	WiggleIterator * source;
} StatData;

static void StatisticSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	StatData * data = (StatData *) wi->data;
	seek(data->source, chrom, start, finish);
	pop(wi);
}

static WiggleIterator * newStatisticIterator(void * data, void (*popFunction)(WiggleIterator *), void (*seekFunction)(WiggleIterator *, const char *, int, int), double default_value, WiggleIterator * source) {
	WiggleIterator * new = newWiggleIterator(data, popFunction, seekFunction, default_value);
	new->append = source;
	return new;
}

//////////////////////////////////////////////////////
// Mean
//////////////////////////////////////////////////////

typedef struct meanData {
	double res;
	double sum;
	double span;
	WiggleIterator * source;
} MeanData;

static void MeanPop(WiggleIterator * wi) {
	MeanData * data = (MeanData *) wi->data;

	if (data->source->done) {
		if (data->span > 0)
			data->res = data->sum / data->span;
		wi->done = true;
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (!isnan(wi->value)) {
		data->sum += (wi->finish - wi->start) * wi->value;
		data->span += (wi->finish - wi->start);
	}
	pop(data->source);
}

static void MeanSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	MeanData * data = (MeanData *) wi->data;
	seek(data->source, chrom, start, finish);
	pop(wi);
}

WiggleIterator * MeanIntegrator(WiggleIterator * wi) {
	MeanData * data = (MeanData *) calloc(1, sizeof(MeanData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->sum = 0;
	data->span = 0;
	data->res = NAN;
	return newStatisticIterator(data, MeanPop, MeanSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// AUC 
//////////////////////////////////////////////////////

static void AUCPop(WiggleIterator * wi) {
	StatData * data = (StatData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (!isnan(wi->value))
		data->res += (wi->finish - wi->start) * wi->value;

	pop(data->source);
}

WiggleIterator * AUCIntegrator(WiggleIterator * wi) {
	StatData * data = (StatData *) calloc(1, sizeof(StatData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = 0;
	return newStatisticIterator(data, AUCPop, StatisticSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Span 
//////////////////////////////////////////////////////

static void SpanPop(WiggleIterator * wi) {
	StatData * data = (StatData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (!isnan(wi->value))
		data->res += (wi->finish - wi->start);

	pop(data->source);
}

WiggleIterator * SpanIntegrator(WiggleIterator * wi) {
	StatData * data = (StatData *) calloc(1, sizeof(StatData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = 0;
	return newStatisticIterator(data, SpanPop, StatisticSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Max 
//////////////////////////////////////////////////////

static void MaxPop(WiggleIterator * wi) {
	StatData * data = (StatData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (!isnan(wi->value) && (isnan(data->res) || wi->value > data->res))
		data->res = wi->value;

	pop(data->source);
}

WiggleIterator * MaxIntegrator(WiggleIterator * wi) {
	StatData * data = (StatData *) calloc(1, sizeof(StatData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = NAN;
	return newStatisticIterator(data, MaxPop, StatisticSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Min 
//////////////////////////////////////////////////////

static void MinPop(WiggleIterator * wi) {
	StatData * data = (StatData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (!isnan(wi->value) && (isnan(data->res) || wi->value < data->res))
		data->res = wi->value;

	pop(data->source);
}

WiggleIterator * MinIntegrator(WiggleIterator * wi) {
	StatData * data = (StatData *) calloc(1, sizeof(StatData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = NAN;
	return newStatisticIterator(data, MinPop, StatisticSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Variance
// Online algorithm copied from 
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass_algorithm
//////////////////////////////////////////////////////

typedef struct varianceData {
	double res;
	WiggleIterator * source;
	double sumWeight;
	double mean;
	double M2;
	double count;
} VarianceData;

static void VariancePop(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		data->res = (data->M2 * data->count) / (data->sumWeight * (data->count - 1));
		return;
	}

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (isnan(wi->value))
		return;
	
	double weight = wi->finish - wi->start;
	double temp = data->sumWeight + weight;
	double delta = wi->value - data->mean;
	double R = delta * weight / temp;
	data->mean += R;
	data->M2 += data->sumWeight * delta * R;
	data->sumWeight = temp;
	data->count++;
	
	pop(data->source);
}

static void VarianceSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	VarianceData * data = (VarianceData *) wi->data;
	seek(data->source, chrom, start, finish);
	pop(wi);
}

WiggleIterator * VarianceIntegrator(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) calloc(1, sizeof(VarianceData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = NAN;
	data->sumWeight = 0;
	data->mean = 0;
	data->M2 = 0;
	data->count = 0;
	return newStatisticIterator(data, VariancePop, VarianceSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Pearson Correlation
//////////////////////////////////////////////////////
// Note: this is an approximate calculation of the Pearson calculation
// which has the benefit of running in a single pass through the data
// The origin of the code can be found at:
// http://en.wikipedia.org/wiki?title=Talk:Correlation

typedef struct pearsonData {
	double res;
	Multiplexer * multi;
	int totalLength;
	double sum_sq_A;
	double sum_sq_B;
	double sum_AB;
	double meanA;
	double meanB;
} PearsonData;

static void PearsonPop(WiggleIterator * wi) {
	PearsonData * data = (PearsonData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		data->res = data->sum_AB / (sqrt(data->sum_sq_A) * sqrt(data->sum_sq_B));
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	wi->value = multi->values[1];

	if (isnan(multi->values[0]) || isnan(multi->values[1]))
		return;

	int halfway, width;
	double sweep, deltaA, deltaB;

	width = (multi->finish - multi->start);
	halfway = data->totalLength + width/2;
	if (halfway == 0)
		halfway = 1;
	sweep = (halfway - 1.0) / halfway;
	data->totalLength += width;
	if (multi->inplay[0])
		deltaA = multi->values[0] - data->meanA;
	else
		deltaA = multi->iters[0]->default_value - data->meanA;
	if (multi->inplay[1])
		deltaB = multi->values[1] - data->meanB;
	else
		deltaB = multi->iters[1]->default_value - data->meanB;
	data->sum_sq_A += deltaA * deltaA * sweep;
	data->sum_sq_B += deltaB * deltaB * sweep;
	data->sum_AB += deltaA * deltaB * sweep;
	data->meanA += deltaA / halfway;
	data->meanB += deltaB / halfway;
	popMultiplexer(multi);
}

WiggleIterator * PearsonIntegrator(WiggleIterator * iterA, WiggleIterator * iterB) {
	PearsonData * data = (PearsonData *) calloc(1, sizeof(PearsonData));
	WiggleIterator * iters[2];
	iters[0] = NonOverlappingWiggleIterator(iterA);
	iters[1] = NonOverlappingWiggleIterator(iterB);
	data->multi = newMultiplexer(iters, 2, false);
	data->totalLength = 0;
	data->sum_sq_A = 0;
	data->sum_sq_B = 0;
	data->sum_AB = 0;
	data->meanA = iterA->value;
	data->meanB = iterB->value;
	data->res = NAN;
	return newStatisticIterator(data, PearsonPop, StatisticSeek, iterB->default_value, iterB);
}

