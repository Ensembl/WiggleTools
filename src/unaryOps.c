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

//////////////////////////////////////////////////////
// Null operator
//////////////////////////////////////////////////////

void NullWiggleIteratorPop(WiggleIterator * wi) {
	return;
}

void NullWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	return;
}

WiggleIterator * NullWiggleIterator() {
	WiggleIterator * new = newWiggleIterator(NULL, &NullWiggleIteratorPop, &NullWiggleIteratorSeek, 0);
	new->done = true;
	return new;
}

//////////////////////////////////////////////////////
// Generic unary operator stuff
//////////////////////////////////////////////////////

typedef struct UnaryWiggleIteratorData_st {
	WiggleIterator * iter;
} UnaryWiggleIteratorData;

void UnaryWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	wi->done = false;
	pop(wi);
}

//////////////////////////////////////////////////////
// Union operator
//////////////////////////////////////////////////////

void UnionWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	int count = 0;

	if (iter->done) {
		wi->done = true;
		return;
	}

	while (!iter->done) {
		if (!count) {
			wi->chrom = iter->chrom;
			wi->start = iter->start;
			wi->finish = iter->finish;
		} else if (wi->chrom == iter->chrom && wi->finish >= iter->start) {
			if (iter->finish > wi->finish)
				wi->finish = iter->finish;
		} else 
			break;
		count++;
		pop(iter);
	} 

	wi->done = (count == 0);
}

WiggleIterator * UnionWiggleIterator(WiggleIterator * i) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &UnionWiggleIteratorPop, &UnaryWiggleIteratorSeek, 0);
}

//////////////////////////////////////////////////////
// Non-overlapping operator
//////////////////////////////////////////////////////

WiggleIterator * NonOverlappingWiggleIterator(WiggleIterator * i) {
	if (i->overlaps)
		return UnionWiggleIterator(i);
	else
		return i;
}

//////////////////////////////////////////////////////
// TestNonOverlapping operator
//////////////////////////////////////////////////////

void TestNonOverlappingWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done) {
		wi->done = true;
	} else {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;

		pop(iter);

		if (wi->chrom == iter->chrom && wi->finish >= iter->start)
			exit(1);
	}
}

WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator * i) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &TestNonOverlappingWiggleIteratorPop, &UnaryWiggleIteratorSeek, i->default_value);
}

//////////////////////////////////////////////////////
// isZero 
//////////////////////////////////////////////////////

void isZeroPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	if (data->iter->done) {
		wi->done = true;
		return;
	}
	wi->chrom = data->iter->chrom;
	wi->start = data->iter->start;
	wi->finish = data->iter->finish;
	wi->value = data->iter->value;
	if (wi->value != 0) {
		exit(1);
	}
	pop(data->iter);
}

WiggleIterator * IsZero(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = wi;
	return newWiggleIterator(data, &isZeroPop, &UnaryWiggleIteratorSeek, wi->default_value);
}

//////////////////////////////////////////////////////
// Default value operator
//////////////////////////////////////////////////////

void DefaultValueWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done) {
		wi->done = true;
	} else {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;
		pop(iter);
	}
}

WiggleIterator * DefaultValueWiggleIterator(WiggleIterator * i, double value) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &DefaultValueWiggleIteratorPop, &UnaryWiggleIteratorSeek, value);
}

//////////////////////////////////////////////////////
// Compression operator
//////////////////////////////////////////////////////

void CompressionWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done) {
		wi->done = true;
	} else {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;
		pop(iter);

		while (!iter->done && strcmp(iter->chrom, wi->chrom) == 0 && iter->start == wi->finish && ((isnan(iter->value) && isnan(wi->value)) || iter->value == wi->value)) {
			wi->finish = iter->finish;
			pop(iter);
		}
	}
}

WiggleIterator * CompressionWiggleIterator(WiggleIterator * i) {
	if (i->overlaps)
		return i;
	else {
		UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
		data->iter = i;
		return newWiggleIterator(data, &CompressionWiggleIteratorPop, &UnaryWiggleIteratorSeek, i->default_value);
	}
}

//////////////////////////////////////////////////////
// Unit operator
//////////////////////////////////////////////////////

void UnitWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!data->iter->done) {
		while (!data->iter->done && (data->iter->value == 0 || isnan(data->iter->value)))
			pop(iter);
		if (data->iter->done) {
			wi->done = true;
			return;
		}
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * UnitWiggleIterator(WiggleIterator * i) {
	if (i->overlaps) {
		return UnionWiggleIterator(i);
	} else {
		UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
		data->iter = i;
		return UnionWiggleIterator(newWiggleIterator(data, &UnitWiggleIteratorPop, &UnaryWiggleIteratorSeek, 0));
	}
}

//////////////////////////////////////////////////////
// High pass filter operator
//////////////////////////////////////////////////////

typedef struct highPassFilterWiggleIteratorData_st {
	WiggleIterator * iter;
	double scalar;
} HighPassFilterWiggleIteratorData;

void HighPassFilterWiggleIteratorPop(WiggleIterator * wi) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!data->iter->done) {
		while (!data->iter->done && (data->iter->value <= data->scalar || isnan(data->iter->value)))
			pop(data->iter);
		if (data->iter->done) {
			wi->done = true;
			return;
		}
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		pop(data->iter);
	} else {
		wi->done = true;
	}
}

void HighPassFilterWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

WiggleIterator * HighPassFilterWiggleIterator(WiggleIterator * i, double s) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) calloc(1, sizeof(HighPassFilterWiggleIteratorData));
	data->iter = i;
	data->scalar = s;
	return UnionWiggleIterator(newWiggleIterator(data, &HighPassFilterWiggleIteratorPop, &HighPassFilterWiggleIteratorSeek, 0));
}

//////////////////////////////////////////////////////
// Overlaps operator
//////////////////////////////////////////////////////

typedef struct overlapWiggleIteratorData_st {
	WiggleIterator * source;
	WiggleIterator * mask;
} OverlapWiggleIteratorData;

void OverlapWiggleIteratorPop(WiggleIterator * wi) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) wi->data;
	WiggleIterator * source = data->source;
	WiggleIterator * mask = data->mask;

	while (!source->done && !mask->done) {
		int chrom_cmp = strcmp(mask->chrom, source->chrom);
		if (chrom_cmp < 0)
			pop(mask);
		else if (chrom_cmp > 0)
			pop(source);
		else if (mask->finish <= source->start)
			pop(mask);
		else if (source->finish <= mask->start)
			pop(source);
		else
			break;
	} 
	
	if (source->done || mask->done)
		wi->done = true;
	else {
		wi->chrom = source->chrom;
		wi->start = source->start;
		wi->finish = source->finish;
		wi->value = source->value;
		pop(source);
	}
}

void OverlapWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) wi->data;
	seek(data->source, chrom, start, finish);
	seek(data->mask, chrom, start, finish);
	pop(wi);
}

WiggleIterator * OverlapWiggleIterator(WiggleIterator * source, WiggleIterator * mask) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) calloc(1, sizeof(OverlapWiggleIteratorData));
	data->source = source;
	data->mask = mask;
	WiggleIterator * wi = newWiggleIterator(data, &OverlapWiggleIteratorPop, &OverlapWiggleIteratorSeek, source->default_value);
	wi->overlaps = source->overlaps;
	return wi;
}

//////////////////////////////////////////////////////
// Scaling operator
//////////////////////////////////////////////////////

typedef struct scaleWiggleIteratorData_st {
	WiggleIterator * iter;
	double scalar;
} ScaleWiggleIteratorData;

void ScaleWiggleIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		if (isnan(iter->value))
			wi->value = NAN;
		else
			wi->value = data->scalar * iter->value;
		pop(data->iter);
	} else {
		wi->done = true;
	}
}

void ScaleWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

WiggleIterator * ScaleWiggleIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->scalar = s;
	float default_value;
	if (isnan(i->default_value))
		default_value = NAN;
	else
		default_value = i->default_value * s;
	return newWiggleIterator(data, &ScaleWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Shifting operator
//////////////////////////////////////////////////////

void ShiftWiggleIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = data->scalar + iter->value;
		pop(data->iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * ShiftWiggleIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->scalar = s;
	float default_value;
	if (isnan(i->default_value))
		default_value = NAN;
	else
		default_value = i->default_value + s;
	return newWiggleIterator(data, &ShiftWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Log operator
//////////////////////////////////////////////////////

const double E = 2.71828128459045;

typedef struct logWiggleIteratorData {
	WiggleIterator * iter;
	double base;
	double baseLog;
} LogWiggleIteratorData;

void LogWiggleIteratorPop(WiggleIterator * wi) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	while (!iter->done && iter->value <= 0)
		pop(iter);

	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		if (isnan(iter->value) || iter->value < 0)
			wi->value = NAN;
		else
			wi->value = log(iter->value) / data->baseLog;
		pop(data->iter);
	} else {
		wi->done = true;
	}
}

void LogWiggleIteratorSeek(WiggleIterator * wi, const char  * chrom, int start, int finish) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

WiggleIterator * NaturalLogWiggleIterator(WiggleIterator * i) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) calloc(1, sizeof(LogWiggleIteratorData));
	data->iter = i;
	data->base = E;
	data->baseLog = 1;
	double default_value;
	if (!isnan(i->default_value) && i->default_value > 0)
		default_value =  log(i->default_value) / data->baseLog;
	else
		default_value = NAN;
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek, default_value);
}

WiggleIterator * LogWiggleIterator(WiggleIterator * i, double s) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) calloc(1, sizeof(LogWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->base = s;
	data->baseLog = log(s);
	double default_value;
	if (!isnan(i->default_value) && i->default_value > 0)
		default_value =  log(i->default_value) / data->baseLog;
	else
		default_value = NAN;
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Exponential operator
//////////////////////////////////////////////////////

typedef struct expWiggleIteratorData_st {
	WiggleIterator * iter;
	double radix;
	double radixLog;
} ExpWiggleIteratorData;

void ExpWiggleIteratorPop(WiggleIterator * wi) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = exp(iter->value * data->radixLog);
		pop(iter);
	} else {
		wi->done = true;
	}
}

void ExpWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

WiggleIterator * ExpWiggleIterator(WiggleIterator * i, double s) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) calloc(1, sizeof(ExpWiggleIteratorData));
	data->iter = i;
	data->radix = s;
	data->radixLog = log(data->radix);
	float default_value;
	if (isnan(i->default_value))
		default_value = NAN;
	else
		default_value = exp(i->default_value * data->radixLog);
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek, default_value);
}

WiggleIterator * NaturalExpWiggleIterator(WiggleIterator * i) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) calloc(1, sizeof(ExpWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->radix = E;
	data->radixLog = 1;
	float default_value;
	if (isnan(i->default_value))
		default_value = NAN;
	else
		default_value = exp(i->default_value * data->radixLog);
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Power operator
//////////////////////////////////////////////////////

static void PowerWiggleIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	// Avoiding divisions by 0
	if (data->scalar < 0)
		while (!iter->done && iter->value == 0)
			pop(iter);

	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		if ((data->scalar < 0 && iter->value <= 0) || isnan(iter->value))
			wi->value = NAN;
		else
			wi->value = pow(iter->value, data->scalar);
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * PowerWiggleIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->scalar = s;
	double default_value;
	if (!isnan(i->default_value) && (i->default_value > 0 || s > 0))
		default_value = pow(i->default_value, s);
	else
		default_value = NAN;
	return newWiggleIterator(data, &PowerWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Abs operator
//////////////////////////////////////////////////////

static void AbsWiggleIteratorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		if (!isnan(iter->value))
			wi->value = abs(iter->value);
		else
			wi->value = NAN;
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * AbsWiggleIterator(WiggleIterator * i) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	double default_value;
	if (!isnan(i->default_value))
		default_value = abs(i->default_value);
	else
		default_value = NAN;
	return newWiggleIterator(data, &AbsWiggleIteratorPop, &UnaryWiggleIteratorSeek, default_value);
}

//////////////////////////////////////////////////////
// Smooth' operator !
//////////////////////////////////////////////////////

typedef struct SmoothWiggleIteratorData_st {
	WiggleIterator * iter;
	double * buffer;
	double sum;
	int first;
	int last;
	int count;
	int width;
} SmoothWiggleIteratorData;
   
static void smoothWiggleIteratorEraseOne(SmoothWiggleIteratorData * data) {
	data->sum -= data->buffer[data->last];
	data->count--;
	if (++data->last >= data->width)
		data->last = 0;
}

static void smoothWiggleIteratorReadOne(char * chrom, int position, SmoothWiggleIteratorData * data) {
	WiggleIterator * iter = data->iter;

	// Let iter run if necessary
	while (!iter->done && iter->chrom == chrom && iter->finish <= position)
		pop(iter);
	
	if (!iter->done && iter->chrom == chrom) {
		// Clear space if necessary
		if (data->count == data->width) 
			smoothWiggleIteratorEraseOne(data);

		// Record iter's value as appropriate
		if (!iter->done && iter->chrom == chrom && iter->start <= position) {
			data->buffer[data->first] = iter->value;
			data->sum += iter->value;
		} else 
			data->buffer[data->first] = 0;
		data->count++;

		// Increment ptr
		if (++data->first >= data->width)
			data->first = 0;

		if (iter->finish == position + 1)
			pop(iter);
	} else 
		smoothWiggleIteratorEraseOne(data);
}

static void SmoothWiggleIteratorPop(WiggleIterator * wi) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	int i;

	if (data->count == 0) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = wi->start + 1;
		for (i = 0; i < data->width / 2 + 1; i++)
			smoothWiggleIteratorReadOne(wi->chrom, wi->start + i, data);
		wi->value = data->sum / data->count;
	} else if (data->count == data->width/2 + 1) {
		if (iter->done) {
			wi->done = true;
			return;
		} else if (iter->chrom != wi->chrom) {
			wi->chrom = iter->chrom;
			wi->start = iter->start;
			wi->finish = wi->start + 1;
			for (i = 0; i < data->width / 2 + 1; i++)
				smoothWiggleIteratorReadOne(wi->chrom, wi->start + i, data);
			wi->value = data->sum / data->count;
		} else {
			wi->start++;
			wi->finish = wi->start + 1;
			if (!iter->done && iter->chrom == wi->chrom)
				smoothWiggleIteratorReadOne(wi->chrom, wi->start + data->width/2, data);
			else
				smoothWiggleIteratorEraseOne(data);
			wi->value = data->sum / data->count;
		}
	} else {
		wi->start++;
		wi->finish = wi->start + 1;
		if (!iter->done && iter->chrom == wi->chrom)
			smoothWiggleIteratorReadOne(wi->chrom, wi->start + data->width/2, data);
		else
			smoothWiggleIteratorEraseOne(data);
		wi->value = data->sum / data->count;
	}
}

void SmoothWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	int i;
	for (i = 0; i < data->width; i++)
		data->buffer[i] = 0;
	data->first = 0;
	data->last = 0;
	data->sum = 0;
	data->count = 0;
	wi->done = false;
	pop(wi);
}

WiggleIterator * SmoothWiggleIterator(WiggleIterator * i, int width) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) calloc(1, sizeof(SmoothWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	data->buffer = (double*) calloc(sizeof(double), width);
	data->width = width;
	return newWiggleIterator(data, &SmoothWiggleIteratorPop, &SmoothWiggleIteratorSeek, i->default_value);
}

//////////////////////////////////////////////////////
// Convenience file reader
//////////////////////////////////////////////////////

WiggleIterator * SmartReader(char * filename) {
	size_t length = strlen(filename);
	if (!strcmp(filename + length - 3, ".bw"))
		return BigWiggleReader(filename);
	else if (!strcmp(filename + length - 3, ".bg"))
		return WiggleReader(filename);
	else if (!strcmp(filename + length - 4, ".wig"))
		return WiggleReader(filename);
	else if (!strcmp(filename + length - 4, ".bed"))
		return BedReader(filename);
	else if (!strcmp(filename + length - 3, ".bb"))
		return BigBedReader(filename);
	else if (!strcmp(filename + length - 4, ".bam"))
		return BamReader(filename);
	else if (!strcmp(filename + length - 4, ".vcf"))
		return VcfReader(filename);
	else if (!strcmp(filename + length - 4, ".bcf"))
		return BcfReader(filename);
	else if (!strcmp(filename, "-"))
		return WiggleReader(filename);
	else {
		fprintf(stderr, "Could not recognize file format from suffix: %s\n", filename);
		exit(1);
	}
}

//////////////////////////////////////////////////////
// Concatenation 
//////////////////////////////////////////////////////

typedef struct CatWiggleIteratorData_st {
	char ** filenames;
	int count;
	int index;
	WiggleIterator * iter;
} CatWiggleIteratorData;

void CatWiggleIteratorPop(WiggleIterator * wi) {
	CatWiggleIteratorData * data = (CatWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;
		pop(iter);
	} else if (data->index < data->count - 1) {
		while (++data->index < data->count) {
			iter = data->iter = SmartReader(data->filenames[data->index]);
			while (!iter->done && (strcmp(wi->chrom, iter->chrom) >= 0 || (strcmp(wi->chrom, iter->chrom) == 0 && wi->finish >= iter->finish)))
				pop(iter);
			if (!iter->done) {
				if (strcmp(wi->chrom, iter->chrom) < 0 || iter->start > wi->finish)
					wi->start = iter->start;
				else
					wi->start = wi->finish;
				wi->chrom = iter->chrom;
				wi->finish = iter->finish;
				wi->value = iter->value;
				pop(iter);
				break;
			}
		}
	} else 
		wi->done = true;
}

void CatWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	puts("Cannot apply seek to a concatenation of files!");
	exit(1);
}

WiggleIterator * CatWiggleIterator(char ** filenames, int count) {
	CatWiggleIteratorData * data = (CatWiggleIteratorData *) calloc(1, sizeof(CatWiggleIteratorData));
	data->count = count;
	data->filenames = filenames;
	data->iter = SmartReader(data->filenames[0]);
	return newWiggleIterator(data, &CatWiggleIteratorPop, &CatWiggleIteratorSeek, 0);
}
