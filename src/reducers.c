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

#include <stdlib.h>
#include <math.h>

#include "multiplexer.h"

typedef struct wiggleReducerData_st {
	Multiplexer * multi;
	bool trim;
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

	while (multi->inplay[data->index]) {
		popMultiplexer(multi);
		if (multi->done) {
			wi->done = true;
			return;
		}
	}

	wi->value = multi->values[data->index];
	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	popMultiplexer(multi);
}

WiggleIterator * SelectReduction(Multiplexer * multi, int index) {
	WiggleSelectData * data = (WiggleSelectData *) calloc(1, sizeof(WiggleSelectData));
	data->multi = multi;
	data->index = index;
	return newWiggleIterator(data, &SelectReductionPop, &WiggleReducerSeek, multi->default_values[index], false);
}

////////////////////////////////////////////////////////
// Fill in 
////////////////////////////////////////////////////////

void FillInReductionPop(WiggleIterator * wi) {
	if (wi->done)
		return;

	WiggleReducerData * data = (WiggleReducerData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		return;
	}

	if (data->trim) {
		while (!multi->inplay[0]) {
			popMultiplexer(multi);
			if (multi->done) {
				wi->done = true;
				return;
			}
		}
	}
	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	if (multi->inplay[1])
		wi->value = multi->values[1];
	else
		wi->value = multi->default_values[1];
	popMultiplexer(multi);
}

WiggleIterator * FillInReduction(Multiplexer * multi, bool trim) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	if (multi->count != 2) {
		printf("The fill in operator can only work on 2 iterators! Got %i\n", multi->count);
		exit(1);
	}
	data->multi = multi;
	data->trim = trim;
	WiggleIterator * res = newWiggleIterator(data, &FillInReductionPop, &WiggleReducerSeek, multi->default_values[1], false);
	return res;
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
		wi->value = multi->values[0];
	else
		wi->value = 0;

	if (isnan(wi->value)) {
		popMultiplexer(multi);
		return;
	}

	for (i = 1; i < multi->count; i++) {
		double value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			break;
		}

		if (wi->value < value)
			wi->value = value;
	}
	popMultiplexer(multi);
}

WiggleIterator * MaxReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double max = data->multi->default_values[0];
	if (!isnan(max)) {
		for (i = 1; i < multi->count; i++) {
			if (isnan(data->multi->default_values[i])) {
				max = NAN;
				break;
			}
			if (data->multi->default_values[i] > max)
				max = data->multi->default_values[i];
		}
	}
	return newWiggleIterator(data, &MaxReductionPop, &WiggleReducerSeek, max, false);
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
		wi->value = multi->values[0];
	else
		wi->value = 0;

	if (isnan(wi->value)) {
		popMultiplexer(multi);
		return;
	}

	for (i = 1; i < multi->count; i++) {
		double value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			break;
		}

		if (wi->value > value)
			wi->value = value;
	}
	popMultiplexer(multi);
}

WiggleIterator * MinReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double min = data->multi->default_values[0];
	if (!isnan(min)) {
		for (i = 1; i < multi->count; i++) {
			if (isnan(data->multi->default_values[i])) {
				min = NAN;
				break;
			}
			if (data->multi->default_values[i] < min)
				min = data->multi->default_values[i];
		}
	}
	return newWiggleIterator(data, &MinReductionPop, &WiggleReducerSeek, min, false);
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
	for (i = 0; i < multi->count; i++) {
		double value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			break;
		}

		wi->value += value;
	}
	popMultiplexer(multi);
}

WiggleIterator * SumReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double sum = 0;
	for (i = 0; i < multi->count; i++) {
		if (isnan(data->multi->default_values[i])) {
			sum = NAN;
			break;
		}
		sum += data->multi->default_values[i];
	}
	return newWiggleIterator(data, &SumReductionPop, &WiggleReducerSeek, sum, false);
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
		double value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			break;
		}

		wi->value *= value;
	}
	popMultiplexer(multi);
}

WiggleIterator * ProductReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double prod = 1;
	for (i = 0; i < multi->count; i++) {
		if (isnan(data->multi->default_values[i])) {
			prod = NAN;
			break;
		}
		prod *= data->multi->default_values[i];
	}
	return newWiggleIterator(data, &ProductReductionPop, &WiggleReducerSeek, prod, false);
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
	for (i = 0; i < multi->count; i++) {
		double value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			break;
		}

		wi->value += value;
	}
	if (!isnan(wi->value))
		wi->value /= multi->count;
	popMultiplexer(multi);
}

WiggleIterator * MeanReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double sum = 0;
	for (i = 0; i < multi->count; i++) {
		if (isnan(data->multi->default_values[i])) {
			sum = NAN;
			break;
		}
		sum += data->multi->default_values[i];
	}
	float default_value;
	if (isnan(sum))
		default_value = NAN;
	else
		default_value = sum/multi->count;
	return newWiggleIterator(data, &MeanReductionPop, &WiggleReducerSeek, default_value, false);
}

////////////////////////////////////////////////////////
// Variance
////////////////////////////////////////////////////////

void VarianceReductionPop(WiggleIterator * wi) {
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

	double mean = 0;
	double count = 0;
	for (i = 0; i < multi->count; i++) {
		float value;

		if (multi->inplay[i])
			value = multi->values[i];
		else 	
			value = multi->default_values[i];

		if (isnan(value)) {
			mean = NAN;
			break;
		}
		mean += value;
		count++;
	}

	if (count < 2 || isnan(mean)) {
		wi->value = NAN;
	} else {
		mean /= count;
		
		wi->value = 0;
		for (i = 0; i < multi->count; i++) {
			if (multi->inplay[i]) {
				double diff = mean - multi->values[i];
				wi->value += diff * diff;
			}
		}
		wi->value /= count;
	}
	popMultiplexer(multi);
}

WiggleIterator * VarianceReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double sum = 0;
	for (i = 0; i < multi->count; i++) {
		if (isnan(data->multi->default_values[i])) {
			sum = NAN;
			break;
		}
		sum += data->multi->default_values[i];
	}
	double default_value;
	if (isnan(sum)) 
		default_value = NAN;
	else {
		double mean = sum / multi->count;
		double error = 0;
		for (i = 0; i < multi->count; i++)
			error += (data->multi->default_values[i] - mean) * (data->multi->default_values[i] - mean);
		default_value = error/multi->count;
	}

	return newWiggleIterator(data, &VarianceReductionPop, &WiggleReducerSeek, default_value, false);
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
	for (i = 0; i < multi->count; i++) {
		float value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];

		if (isnan(value)) {
			mean = NAN;
			break;
		}

		mean += value;
	}

	if (isnan(mean))
		wi->value = NAN;
	else {
		mean /= multi->count;
		
		wi->value = 0;
		for (i = 0; i < multi->count; i++) {
			if (multi->inplay[i])
				diff = mean - multi->values[i];
			else
				diff = mean - multi->default_values[i];
			wi->value += diff * diff;
		}
		wi->value /= multi->count;
		wi->value = sqrt(wi->value);
	}
	
	popMultiplexer(multi);
}

WiggleIterator * StdDevReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;
	int i;
	double sum = 0;
	for (i = 0; i < multi->count; i++) {
		if (isnan(data->multi->default_values[i])) {
			sum = NAN;
			break;
		}
		sum += data->multi->default_values[i];
	}
	double default_value;

	if (isnan(sum)) 
		default_value = NAN;
	else {
		double mean = sum / multi->count;
		double error = 0;
		for (i = 0; i < multi->count; i++)
			error += (data->multi->default_values[i] - mean) * (data->multi->default_values[i] - mean);
		default_value = sqrt(error/multi->count);
	}

	return newWiggleIterator(data, &StdDevReductionPop, &WiggleReducerSeek, default_value, false);
}

////////////////////////////////////////////////////////
// Shannon entropy
////////////////////////////////////////////////////////

void EntropyReductionPop(WiggleIterator * wi) {
	int i;
	float count = 0;

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

	for (i = 0; i < multi->count; i++) {
		double value;
		if (multi->inplay[i]) 
			value = multi->values[i]; 
		else 
			value = multi->default_values[i];

		if (isnan(value)) {
			wi->value = NAN;
			popMultiplexer(multi);
			return;
		} else if (value > 0)
			count++;
	}

	if (count == 0 || count == multi->count)
		wi->value = 0;
	else {
		double p = count / multi->count;
		wi->value = - p * log(p) - (1-p) * log(1 - p);
	}
	
	popMultiplexer(multi);
}

WiggleIterator * EntropyReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;

	int i;
	int count = 0;
	for (i = 0; i < multi->count; i++) {
		if (isnan(multi->default_values[i])) {
			count = -1;
			break;
		}
		if (multi->default_values[i] != 0)
			count++;
	}
	double default_value;
	if (count == -1)
		default_value = NAN;
	else {
		double p = count / multi->count;
		if (p)
			default_value = - p * log(p) - (1-p) * log(1 - p);
		else
			default_value = 0;
	}

	return newWiggleIterator(data, &StdDevReductionPop, &WiggleReducerSeek, default_value, false);
}

////////////////////////////////////////////////////////
// CV
////////////////////////////////////////////////////////

void CVReductionPop(WiggleIterator * wi) {
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
	for (i = 0; i < multi->count; i++) {
		float value;
		if (multi->inplay[i])
			value = multi->values[i];
		else
			value = multi->default_values[i];
		mean += value;

		if (isnan(value)) {
			wi->value = NAN;
			popMultiplexer(multi);
			return;
		}
	}
	mean /= multi->count;

	if (mean == 0) {
		wi->value = NAN;
		popMultiplexer(multi);
		return;
	}
	
	wi->value = 0;
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i])
			diff = mean - multi->values[i];
		else
			diff = mean - multi->default_values[i];
		wi->value += diff * diff;
	}
	wi->value /= multi->count;
	wi->value = sqrt(wi->value) / mean;
	
	popMultiplexer(multi);
}

WiggleIterator * CVReduction(Multiplexer * multi) {
	WiggleReducerData * data = (WiggleReducerData *) calloc(1, sizeof(WiggleReducerData));
	data->multi = multi;

	int i;
	double mean = 0;
	for (i = 0; i < multi->count; i++) {
		if ( isnan(multi->default_values[i])) {
			mean = NAN;
			break;
		} 
		mean += multi->default_values[i];
	}

	float default_value;
	if (!isnan(mean)) {
		mean /= multi->count;
		double error = 0;
		for (i = 0; i < multi->count; i++)
			error += (mean - multi->default_values[i]) * (mean - multi->default_values[i]);
		default_value = sqrt(error/multi->count)/mean;
	} else
		default_value = NAN;
	return newWiggleIterator(data, &CVReductionPop, &WiggleReducerSeek, default_value, false);
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
	pop(iter);
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
			data->vals[i] = multi->default_values[i];

		if (isnan(data->vals[i])) {
			wi->value = NAN;
			popMultiplexer(multi);
			return;
		}
	}
	qsort(data->vals, multi->count, sizeof(double), &compDoubles);
	wi->value = data->vals[multi->count / 2];
	
	popMultiplexer(multi);
}

WiggleIterator * MedianReduction(Multiplexer * multi) {
	MedianWiggleReducerData * data = (MedianWiggleReducerData *) calloc(1, sizeof(MedianWiggleReducerData));
	data->multi = multi;
	data->vals = (double *) calloc(data->multi->count, sizeof(double));
	int i;
	float default_value = 0;
	for (i = 0; i < multi->count; i++) {
		data->vals[i] = multi->default_values[i];
		if (isnan(data->vals[i])) {
			default_value = NAN;
			break;
		}
	}

	if (!isnan(default_value)) {
		qsort(data->vals, multi->count, sizeof(double), &compDoubles);
		default_value = data->vals[multi->count/2];
	}
	return newWiggleIterator(data, &MedianReductionPop, &MedianWiggleReducerSeek, default_value, false);
}
