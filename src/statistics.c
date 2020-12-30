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

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Local header
#include "wiggleIterator.h"
#include "multiplexer.h"
#include "multiSet.h"

#define PI (3.141592653589793)

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
	WiggleIterator * new = newWiggleIterator(data, popFunction, seekFunction, default_value, false);
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
// See Technical supplement for explanations
//////////////////////////////////////////////////////

typedef struct varianceData {
	double res;
	WiggleIterator * source;
	double T;
	double sum;
	long count;
} VarianceData;

static void VarianceCorePop(WiggleIterator * wi, VarianceData * data) {
	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	if (isnan(wi->value))
		return;

	int length = wi->finish - wi->start;
	if (data->count) {
		double old_mean = data->sum / data->count;
		double new_mean = data->sum / (data->count + length);
		double delta_T = old_mean * new_mean - new_mean * 2 * wi->value + ((double) data->count / (data->count + length)) * wi->value * wi->value;
		data->T += delta_T * length;
	}
	data->count += length;
	data->sum += length * wi->value;
	
	pop(data->source);
}

static void VariancePop(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		data->res = data->T / (data->count - 1);
		return;
	}

	VarianceCorePop(wi, data);
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
	return newStatisticIterator(data, VariancePop, VarianceSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Standard deviation 
//////////////////////////////////////////////////////

static void StandardDeviationPop(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		data->res = data->T / (data->count - 1);
		data->res = sqrt(data->res);
		return;
	}

	VarianceCorePop(wi, data);
}

WiggleIterator * StandardDeviationIntegrator(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) calloc(1, sizeof(VarianceData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = NAN;
	return newStatisticIterator(data, StandardDeviationPop, VarianceSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Coefficient of Variation
//////////////////////////////////////////////////////

static void CoefficientOfVariationPop(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		data->res = data->T / (data->count - 1);
		data->res = sqrt(data->res);
		data->res /= (data->sum / data->count);
		return;
	}

	VarianceCorePop(wi, data);
}

WiggleIterator * CoefficientOfVariationIntegrator(WiggleIterator * wi) {
	VarianceData * data = (VarianceData *) calloc(1, sizeof(VarianceData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->res = NAN;
	return newStatisticIterator(data, CoefficientOfVariationPop, VarianceSeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Energy
//////////////////////////////////////////////////////
// Computes the square norm of the Fourier transform at
// a given wavelength
//
// Note: this simply computes the value of the Fourier
// transform at a single wavelength. Computing the entire
// Fourier transform of an input signal would require 
//
// a) FFT implementation (e.g. KissFFT:
// https://sourceforge.net/projects/kissfft/)
// 
// b) Some very clever memomry management to store a massive
// DFT matrix
//

typedef struct energyData_st {
	double res;
	double real;
	double im;
	int wavelength;
	int chrom_offset;
	WiggleIterator * source;
} EnergyData;

static void EnergySeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	EnergyData * data = (EnergyData *) wi->data;
	seek(data->source, chrom, start, finish);
	pop(wi);
}

static void EnergyPop(WiggleIterator * wi) {
	EnergyData * data = (EnergyData *) wi->data;

	if (data->source->done) {
		wi->done = true;
		data->res = data->real * data->real + data->im * data->im;
		return;
	}

	if (wi->chrom != data->source->chrom)
		data->chrom_offset += wi->finish;

	wi->chrom = data->source->chrom;
	wi->start = data->source->start;
	wi->finish = data->source->finish;
	wi->value = data->source->value;

	int position;
	for(position = data->chrom_offset + data->source->start; position < data->chrom_offset + data->source->finish; position++) {
		data->real += cos(- position * 2 * PI / data->wavelength) * data->source->value;
		data->im += sin(- position * 2 * PI / data->wavelength) * data->source->value;
	}
	pop(data->source);
}

WiggleIterator * EnergyIntegrator(WiggleIterator * wi, int wavelength) {
	EnergyData * data = (EnergyData *) calloc(1, sizeof(EnergyData));
	data->source = NonOverlappingWiggleIterator(wi);
	data->wavelength = wavelength;
	data->res = NAN;
	return newStatisticIterator(data, EnergyPop, EnergySeek, wi->default_value, wi);
}

//////////////////////////////////////////////////////
// Pearson Correlation
//////////////////////////////////////////////////////

typedef struct pearsonData {
	double res;
	Multiplexer * multi;
	int count;
	double sum_X;
	double sum_Y;
	double T_XX;
	double T_XY;
	double T_YY;
} PearsonData;

static void PearsonSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	PearsonData * data = (PearsonData *) wi->data;
	seekMultiplexer(data->multi, chrom, start, finish);
	pop(wi);
}

static void PearsonPop(WiggleIterator * wi) {
	PearsonData * data = (PearsonData *) wi->data;
	Multiplexer * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		if (data->T_XX * data->T_YY)
			data->res = data->T_XY / sqrt(data->T_XX * data->T_YY);
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;
	wi->value = multi->values[1];

	double X, Y;

	if (multi->inplay[0])
		X = multi->values[0];
	else
		X = multi->iters[0]->default_value;

	if (multi->inplay[1])
		Y = multi->values[1];
	else
		Y = multi->iters[1]->default_value;

	int length = (multi->finish - multi->start);
	if (data->count) {
		double old_mean_X = data->sum_X / data->count;
		double new_mean_X = data->sum_X / (data->count + length);
		double old_mean_Y = data->sum_Y / data->count;
		double new_mean_Y = data->sum_Y / (data->count + length);
		double scaling_ratio = (double) data->count / (data->count + length);

		data->T_XY += (new_mean_X * old_mean_Y + scaling_ratio * X * Y - new_mean_X * Y - new_mean_Y * X) * length;
		data->T_XX += (new_mean_X * (old_mean_X - 2 * X) + scaling_ratio * X * X) * length;
		data->T_YY += (new_mean_Y * (old_mean_Y - 2 * Y) + scaling_ratio * Y * Y) * length;
	}
	data->count += length;
	data->sum_X += X * length;
	data->sum_Y += Y * length;
	popMultiplexer(multi);
}

WiggleIterator * PearsonIntegrator(Multiplexer * multi) {
	PearsonData * data = (PearsonData *) calloc(1, sizeof(PearsonData));
	data->multi = multi;
	data->res = NAN;
	return newStatisticIterator(data, PearsonPop, PearsonSeek, data->multi->iters[1]->default_value, data->multi->iters[1]);
}

////////////////////////////////////////////////////////
// N-dimensional Pearson 
////////////////////////////////////////////////////////

typedef struct ndpearsonData_st {
	double res;
	WiggleIterator * source;
	Multiset * multi;
	int rank;
	int count;
	double * sum_X;
	double * sum_Y;
	double T_XX;
	double T_XY;
	double T_YY;
} NDPearsonData;

void NDPearsonSeek(WiggleIterator * iter, const char * chrom, int start, int finish) {
	NDPearsonData * data = (NDPearsonData* ) iter->data;
	seekMultiset(data->multi, chrom, start, finish);
	pop(iter);
}

void NDPearsonPop(WiggleIterator * wi) {
	if (wi->done)
		return;

	NDPearsonData * data = (NDPearsonData *) wi->data;
	Multiset * multi = data->multi;

	if (multi->done) {
		wi->done = true;
		if (data->T_XX * data->T_YY)
			data->res = data->T_XY / sqrt(data->T_XX * data->T_YY);
		return;
	}

	wi->chrom = multi->chrom;
	wi->start = multi->start;
	wi->finish = multi->finish;

	int length = (multi->finish - multi->start);
	double scaling_ratio = (double) data->count / (data->count + length);
	int dim;
	for (dim = 0; dim < data->rank; dim++) {
		double Xi, Yi;

		if (multi->inplay[0])
			Xi = multi->values[0][dim];
		else
			Xi = multi->multis[0]->iters[dim]->default_value;

		if (multi->inplay[1])
			Yi = multi->values[1][dim];
		else
			Yi = multi->multis[1]->iters[dim]->default_value;

		if (data->count) {
			double old_mean_Xi = data->sum_X[dim] / data->count;
			double new_mean_Xi = data->sum_X[dim] / (data->count + length);
			double old_mean_Yi = data->sum_Y[dim] / data->count;
			double new_mean_Yi = data->sum_Y[dim] / (data->count + length);

			data->T_XY += (new_mean_Xi * old_mean_Yi + scaling_ratio * Xi * Yi - new_mean_Xi * Yi - new_mean_Yi * Xi) * length;
			data->T_XX += (new_mean_Xi * (old_mean_Xi - 2 * Xi) + scaling_ratio * Xi * Xi) * length;
			data->T_YY += (new_mean_Yi * (old_mean_Yi - 2 * Yi) + scaling_ratio * Yi * Yi) * length;
		}
		data->sum_X[dim] += Xi * length;
		data->sum_Y[dim] += Yi * length;
	}

	data->count += length;

	// Update inputs
	popMultiset(multi);
}

WiggleIterator * NDPearsonIntegrator(Multiset * multi) {
	if (multi->count != 2 || multi->multis[0]->count != multi->multis[1]->count) {
		fprintf(stderr, "Incorrect number of input tracks to N-dimensional Pearson correlation!\n");
		exit(1);
	}
	NDPearsonData * data = (NDPearsonData *) calloc(1, sizeof(NDPearsonData));
	data->multi = multi;
	data->rank = multi->count;
	data->sum_X = calloc(data->rank, sizeof(double));
	data->sum_Y = calloc(data->rank, sizeof(double));
	data->res = NAN;
	return newStatisticIterator(data, NDPearsonPop, NDPearsonSeek, 0, multi->multis[0]->iters[0]);
}

//////////////////////////////////////////////////////
// Print statistics operator
//////////////////////////////////////////////////////

typedef struct printStatisticsData_st {
	WiggleIterator * iter;
	FILE * file;
} PrintStatisticsData;

static void PrintStatisticsWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	PrintStatisticsData * data = (PrintStatisticsData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

static void PrintStatisticsWiggleIteratorPop(WiggleIterator * wi) {
	PrintStatisticsData * data = (PrintStatisticsData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done) {
		wi->done = true;
		if (iter->append)
			fprintf(data->file, "%f", *((double *) iter->data));
		for (iter = iter->append; iter && iter->append; iter = iter->append)
			fprintf(data->file, "\t%f", *((double *) iter->data));
		fprintf(data->file, "\n");
	} else {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;
		pop(iter);
	}
}

WiggleIterator * PrintStatisticsWiggleIterator(WiggleIterator * i, FILE * file) {
	PrintStatisticsData * data = (PrintStatisticsData *) calloc(1, sizeof(PrintStatisticsData));
	data->iter = i;
	data->file = file;
	return newWiggleIterator(data, &PrintStatisticsWiggleIteratorPop, &PrintStatisticsWiggleIteratorSeek, i->default_value, false);
}
