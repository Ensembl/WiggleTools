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
// Profile summaries
//////////////////////////////////////////////////////

static void updateProfile(WiggleIterator * wig, double compression, double * profile, int profile_width, bool stranded) {
	int start, finish, pos;

	if (!stranded || wig->strand > 0) {
		start = (int) round(wig->start * compression);
		finish = (int) round(wig->finish * compression); 
	} else if (wig->strand < 0) {
		start = (int) round(profile_width - 1 - (wig->finish * compression));
		finish = (int) round(profile_width - 1 - (wig->start * compression)); 
	} else {
		fprintf(stderr, "Cannot provide stranded profile on non-stranded regions\n");
		exit(1);
	}

	if (start < 0)
		start = 0;
	if (finish <= start)
		finish = start + 1;
	if (finish > profile_width)
		finish = profile_width;

	for (pos = start; pos < finish; pos++)
		profile[pos] += wig->value / compression;
}

void regionProfile(WiggleIterator * wig, double * profile, int profile_width, int region_width, bool stranded) {
	double compression = profile_width / (double) region_width;
	int pos;

	for (pos = 0; pos < profile_width; pos++)
		profile[pos] = 0;

	for (; !wig->done; pop(wig))
		updateProfile(wig, compression, profile, profile_width, stranded);
}

void addProfile(double * dest, double * source, int width) {
	int i;

	for (i=0; i<width; i++) 
		dest[i] += source[i];
}

//////////////////////////////////////////////////////
// Histograms
//////////////////////////////////////////////////////

struct histogram_st {
	int width;
	double * values;
	double min;
	double max;
};

static void reassignColumnRight(Histogram * hist, double ratio, int column) {
	double start = hist->width - (hist->width - column) * ratio;
	int int_start = (int) start;
	double end = hist->width - (hist->width - column - 1) * ratio;
	int int_end = (int) end;
	double value = hist->values[column];
	if (int_start == int_end) {
		hist->values[int_start] += value;
	} else {
		double split =  (end - int_end) / ratio;
		hist->values[int_end] += value * split;
		hist->values[int_start] += value * (1 - split); 
	}
	hist->values[column] -= value;
}

static void lowerMin(Histogram * hist, double value) {
	if (hist->max != hist->min) {
		double ratio = (hist->max - hist->min) / (hist->max - value);
		int column;
		for (column = 0; column < hist->width - 1 ; column++) {
			reassignColumnRight(hist, ratio, column);
		}
	} else {
		hist->values[hist->width - 1] = hist->values[0];
		hist->values[0] = 0;
	}
	hist->min = value;
}

static void reassignColumnLeft(Histogram * hist, double ratio, int column) {
	double start = column * ratio;
	int int_start = (int) start;
	double end = (column + 1) * ratio;
	int int_end = (int) end;
	double value = hist->values[column];
	if (int_start == int_end) {
		hist->values[int_start] += value;
	} else {
		double split =  (end - int_end) / ratio;
		hist->values[int_end] += value * split;
		hist->values[int_start] += value * (1 - split); 
	}
	hist->values[column] -= value;
}

static void raiseMax(Histogram * hist, double value) {
	if (hist->max != hist->min) {
		double ratio = (hist->max - hist->min) / (value - hist->min);
		int column;
		for (column = 1; column < hist->width; column++) {
			reassignColumnLeft(hist, ratio, column);
		}
	}
	hist->max = value;
}

static void insertIntoHistogram(Histogram * hist, WiggleIterator * wig) {
	int column = (int) ((wig->value - hist->min) * (hist->width - 1) / (hist->max - hist->min));
	hist->values[column] += wig->finish - wig->start;
}

static void updateHistogram(Histogram * hist, WiggleIterator * wig) {
	if (wig->value > hist->max)
		raiseMax(hist, wig->value);
	else if (wig->value < hist->min)
		lowerMin(hist, wig->value);

	if (hist->min != hist->max)
		insertIntoHistogram(hist, wig);
	else 
		hist->values[0] += wig->finish - wig->start;
}

Histogram * histogram(WiggleIterator * wig, int width) {
	Histogram * hist = calloc(1, sizeof(Histogram));
	hist->width = width;
	hist->values = calloc(width, sizeof(double));
	hist->min = hist->max = wig->value;
	hist->values[0] = wig->finish - wig->start;
	pop(wig);

	for (; !wig->done; pop(wig))
		updateHistogram(hist, wig);

	return hist;
}

void normalize_histogram(Histogram * hist) {
	double sum = 0;
	int column;

	for (column = 0; column < hist->width; column++)
		sum += hist->values[column];

	for (column = 0; column < hist->width; column++)
		hist->values[column] /= sum;
}

void print_histogram(Histogram * hist, FILE * file) {
	double step = (hist->max - hist->min) / hist->width;
	double position = hist->min;
	int column;

	for (column = 0; column < hist->width; column++) {
		fprintf(file, "%f\t%f\n", position, hist->values[column]);
		position += step;
	}
	fprintf(file, "%f\t.\n", position);
}
