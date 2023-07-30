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
#include "fib.h"

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
	WiggleIterator * new = newWiggleIterator(NULL, &NullWiggleIteratorPop, &NullWiggleIteratorSeek, 0, false);
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
			wi->value = iter->value;
		} else if (wi->chrom == iter->chrom && wi->finish > iter->start) {
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
	return newWiggleIterator(data, &UnionWiggleIteratorPop, &UnaryWiggleIteratorSeek, i->default_value, false);
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
	return newWiggleIterator(data, &TestNonOverlappingWiggleIteratorPop, &UnaryWiggleIteratorSeek, i->default_value, false);
}

//////////////////////////////////////////////////////
// floor
//////////////////////////////////////////////////////

void floorPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	if (data->iter->done) {
		wi->done = true;
		return;
	}
	wi->chrom = data->iter->chrom;
	wi->start = data->iter->start;
	wi->finish = data->iter->finish;
    wi->value = floor(data->iter->value);
	pop(data->iter);
}

WiggleIterator * Floor(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = wi;
	return newWiggleIterator(data, &floorPop, &UnaryWiggleIteratorSeek, floor(wi->default_value), wi->overlaps);
}

//////////////////////////////////////////////////////
// toInt
//////////////////////////////////////////////////////

void toIntPop(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) wi->data;
	if (data->iter->done) {
		wi->done = true;
		return;
	}
	wi->chrom = data->iter->chrom;
	wi->start = data->iter->start;
	wi->finish = data->iter->finish;
	wi->value = (int)(data->iter->value);
	pop(data->iter);
}

WiggleIterator * ToInt(WiggleIterator * wi) {
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = wi;
	return newWiggleIterator(data, &toIntPop, &UnaryWiggleIteratorSeek, (int)(wi->default_value), wi->overlaps);
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
	return newWiggleIterator(data, &isZeroPop, &UnaryWiggleIteratorSeek, wi->default_value, wi->overlaps);
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
	return newWiggleIterator(data, &DefaultValueWiggleIteratorPop, &UnaryWiggleIteratorSeek, value, i->overlaps);
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

		while (!iter->done && strcmp(iter->chrom, wi->chrom) == 0 && iter->start == wi->finish && ((isnan(iter->value) && isnan(wi->value)) || (fabs(iter->value - wi->value) < 0.000001))) {
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
		data->iter = NonOverlappingWiggleIterator(i);
		return newWiggleIterator(data, &CompressionWiggleIteratorPop, &UnaryWiggleIteratorSeek, i->default_value, false);
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
	UnaryWiggleIteratorData * data = (UnaryWiggleIteratorData *) calloc(1, sizeof(UnaryWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	return newWiggleIterator(data, &UnitWiggleIteratorPop, &UnaryWiggleIteratorSeek, 0, false);
}

//////////////////////////////////////////////////////
// Coverage operator
//////////////////////////////////////////////////////

typedef struct CoverageWiggleIteratorData_st {
	WiggleIterator * iter;
	FibHeap * heap;
} CoverageWiggleIteratorData;

void CoverageWiggleIteratorPop(WiggleIterator * wi) {
	CoverageWiggleIteratorData * data = (CoverageWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		if (wi->chrom == NULL)
			wi->value = 0;

		while (fh_notempty(data->heap) && fh_min(data->heap) == wi->finish) {
			wi->value--;
			fh_extractmin(data->heap);
		}

		if (wi->value < 0) {
			fprintf(stderr, "Negative coverage???\n");
			exit(1);
		}

		if (wi->value) {
			wi->start = wi->finish;
		} else {
			wi->chrom = iter->chrom;
			wi->start = iter->start;
		}

		while (!iter->done && !strcmp(iter->chrom, wi->chrom) && iter->start == wi->start) {
			fh_insert(data->heap, iter->finish, 0);
			pop(iter);
			wi->value++;
		}

		if (!fh_notempty(data->heap) || (!strcmp(iter->chrom, wi->chrom) && iter->start < fh_min(data->heap)))
			wi->finish = iter->start;
		else
			wi->finish = fh_min(data->heap);
	} else if (fh_notempty(data->heap)) {
		wi->start = wi->finish;
		while (fh_notempty(data->heap) && fh_min(data->heap) == wi->start) {
			wi->value--;
			fh_extractmin(data->heap);
		}

		if (wi->value)
			wi->finish = fh_min(data->heap);
		else {
			wi->done = true;
			fh_deleteheap(data->heap);
			data->heap = NULL;
		}
	} else {
		wi->done = true;
		fh_deleteheap(data->heap);
		data->heap = NULL;
	}
}

void CoverageWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	CoverageWiggleIteratorData * data = (CoverageWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	wi->value = 0;
	fh_deleteheap(data->heap);
	data->heap = fh_makeheap();
	pop(wi);
}

WiggleIterator * CoverageWiggleIterator(WiggleIterator * i) {
	if (i->overlaps) {
		CoverageWiggleIteratorData * data = (CoverageWiggleIteratorData *) calloc(1, sizeof(CoverageWiggleIteratorData));
		data->iter = i;
		data->heap = fh_makeheap();
		return newWiggleIterator(data, &CoverageWiggleIteratorPop, &CoverageWiggleIteratorSeek, 0, false);
	} else
		return i;
}

//////////////////////////////////////////////////////
// High pass filter operator
//////////////////////////////////////////////////////

typedef struct highPassFilterWiggleIteratorData_st {
	WiggleIterator * iter;
	double scalar;
	bool equal;
} HighPassFilterWiggleIteratorData;

void HighPassFilterWiggleIteratorPop(WiggleIterator * wi) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	while (!data->iter->done 
			&& 
			( 
			 (data->equal 
			  && data->iter->value < data->scalar)
			 || (!data->equal
				 && data->iter->value <= data->scalar)
			 || isnan(data->iter->value)
			 )
			)
		pop(data->iter);

	if (data->iter->done) {
		wi->done = true;
		return;
	}

	wi->chrom = iter->chrom;
	wi->start = iter->start;
	wi->finish = iter->finish;
	pop(data->iter);
}

void HighPassFilterWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	pop(wi);
}

WiggleIterator * HighPassFilterWiggleIterator(WiggleIterator * i, double s, bool equal) {
	HighPassFilterWiggleIteratorData * data = (HighPassFilterWiggleIteratorData *) calloc(1, sizeof(HighPassFilterWiggleIteratorData));
	data->iter = i;
	data->scalar = s;
	data->equal = equal;
	return newWiggleIterator(data, &HighPassFilterWiggleIteratorPop, &HighPassFilterWiggleIteratorSeek, 0, i->overlaps);
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
	return newWiggleIterator(data, &OverlapWiggleIteratorPop, &OverlapWiggleIteratorSeek, source->default_value, source->overlaps);
}

//////////////////////////////////////////////////////
// Trim operator
//////////////////////////////////////////////////////

void TrimWiggleIteratorPop(WiggleIterator * wi) {
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
		wi->start = source->start > mask->start? source->start: mask->start;
		wi->finish = source->finish < mask->finish? source->finish: mask->finish;
		wi->value = source->value;
		if (source->finish <= mask->finish)
			pop(source);
		else
			pop(mask);
	}
}

WiggleIterator * TrimWiggleIterator(WiggleIterator * source, WiggleIterator * mask) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) calloc(1, sizeof(OverlapWiggleIteratorData));
	data->source = source;
	data->mask = NonOverlappingWiggleIterator(mask);
	return newWiggleIterator(data, &TrimWiggleIteratorPop, &OverlapWiggleIteratorSeek, source->default_value, source->overlaps);
}

//////////////////////////////////////////////////////
// No-Overlaps operator
//////////////////////////////////////////////////////

void NoverlapWiggleIteratorPop(WiggleIterator * wi) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) wi->data;
	WiggleIterator * source = data->source;
	WiggleIterator * mask = data->mask;

	while (!source->done && !mask->done) {
		int chrom_cmp = strcmp(mask->chrom, source->chrom);
		if (chrom_cmp < 0)
			pop(mask);
		else if (chrom_cmp > 0)
			break;
		else if (mask->finish <= source->start)
			pop(mask);
		else if (source->finish <= mask->start)
			break;
		else
			pop(source);
	}

	if (source->done)
		wi->done = true;
	else {
		wi->chrom = source->chrom;
		wi->start = source->start;
		wi->finish = source->finish;
		wi->value = source->value;
		pop(source);
	}
}

WiggleIterator * NoverlapWiggleIterator(WiggleIterator * source, WiggleIterator * mask) {
	OverlapWiggleIteratorData * data = (OverlapWiggleIteratorData *) calloc(1, sizeof(OverlapWiggleIteratorData));
	data->source = source;
	data->mask = mask;
	return newWiggleIterator(data, &NoverlapWiggleIteratorPop, &OverlapWiggleIteratorSeek, source->default_value, source->overlaps);
}

//////////////////////////////////////////////////////
// Nearest operator
//////////////////////////////////////////////////////

typedef struct nearestWiggleIteratorData_st {
	WiggleIterator * source;
	WiggleIterator * mask;
	char * prev_chrom;
	int prev_start;
	int prev_finish;
} NearestWiggleIteratorData;

void NearestWiggleIteratorPop(WiggleIterator * wi) {
	NearestWiggleIteratorData * data = (NearestWiggleIteratorData *) wi->data;
	WiggleIterator * source = data->source;
	WiggleIterator * mask = data->mask;

	if (source->done) {
		wi->done = true;
		return;
	}

	while (!mask->done) {
		int chrom_cmp = strcmp(mask->chrom, source->chrom);
		if (chrom_cmp < 0)
			pop(mask);
		else if (chrom_cmp > 0)
			break;
		else if (mask->start <= source->start) {
			data->prev_chrom = mask->chrom;
			data->prev_start = mask->start;
			data->prev_finish = mask->finish;
			pop(mask);
		} else
			break;
	}

	wi->chrom = source->chrom;
	wi->start = source->start;
	wi->finish = source->finish;

	bool set = false;
	if (data->prev_chrom && strcmp(data->prev_chrom, source->chrom) == 0) {
		wi->value = wi->start - data->prev_finish + 1;
		set = true;
	}

	if (!mask->done && strcmp(mask->chrom, source->chrom) == 0 && (!set || wi->value > mask->start - wi->finish + 1)) {
		wi->value = mask->start - wi->finish + 1;
		set = true;
	}

	if (!set)
		wi->value = NAN;
	else if (wi->value < 0)
		wi->value = 0;

	pop(source);
}

void NearestWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	NearestWiggleIteratorData * data = (NearestWiggleIteratorData *) wi->data;
	data->prev_chrom = NULL;
	seek(data->source, chrom, start, finish);
	seek(data->mask, chrom, start, finish);
	pop(wi);
}

WiggleIterator * NearestWiggleIterator(WiggleIterator * source, WiggleIterator * mask) {
	NearestWiggleIteratorData * data = (NearestWiggleIteratorData *) calloc(1, sizeof(NearestWiggleIteratorData));
	data->source = source;
	data->mask = mask;
	return newWiggleIterator(data, &NearestWiggleIteratorPop, &NearestWiggleIteratorSeek, source->default_value, source->overlaps);
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
	return newWiggleIterator(data, &ScaleWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value, i->overlaps);
}

//////////////////////////////////////////////////////
// ShiftPos operator
//////////////////////////////////////////////////////

void ShiftPosIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		while (data->scalar >= iter->finish) {
			pop(iter);
		}
		wi->chrom = iter->chrom;
		wi->start = (data->scalar >= iter->start ? 1 : iter->start - data->scalar);
		wi->finish = iter->finish - data->scalar;
		wi->value = iter->value;
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * ShiftPosIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = NonOverlappingWiggleIterator(i);
	if (s < 0) {
		fprintf(stderr, "Cannot provide a negative value. All coordinates are shifted downwards by default.\n");
		exit(1);
	}
	data->scalar = s;
	float default_value = i->default_value;
	return newWiggleIterator(data, &ShiftPosIteratorPop, &ScaleWiggleIteratorSeek, default_value, i->overlaps);
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
	return newWiggleIterator(data, &ShiftWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value, i->overlaps);
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
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek, default_value, i->overlaps);
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
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek, default_value, i->overlaps);
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
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek, default_value, i->overlaps);
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
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek, default_value, i->overlaps);
}

//////////////////////////////////////////////////////
// Power operator
//////////////////////////////////////////////////////

static void PowerWiggleIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

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
	return newWiggleIterator(data, &PowerWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value, i->overlaps);
}

//////////////////////////////////////////////////////
// Extend operator
//////////////////////////////////////////////////////

static void ExtendWiggleIteratorPop(WiggleIterator * wi) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start - data->scalar;
		wi->finish = iter->finish + data->scalar;
		wi->value = iter->value;
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * ExtendWiggleIterator(WiggleIterator * i, int s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = i;
	data->scalar = s;
	double default_value = i->default_value;
	return newWiggleIterator(data, &ExtendWiggleIteratorPop, &ScaleWiggleIteratorSeek, default_value, i->overlaps);
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
			wi->value = fabs(iter->value);
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
		default_value = fabs(i->default_value);
	else
		default_value = NAN;
	return newWiggleIterator(data, &AbsWiggleIteratorPop, &UnaryWiggleIteratorSeek, default_value, i->overlaps);
}

//////////////////////////////////////////////////////
// Binning operator
//////////////////////////////////////////////////////

typedef struct BinningWiggleIteratorData_st {
	WiggleIterator * iter;
	int width;
} BinningWiggleIteratorData;

void BinningWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BinningWiggleIteratorData * data = (BinningWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	wi->chrom = NULL;
	pop(wi);
}

static void BinningWiggleIteratorPop(WiggleIterator * wi) {
	BinningWiggleIteratorData * data = (BinningWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	int width = data->width;

	if (iter->done) {
		wi->done = true;
		return;
	}

	// Set new boundaries
	// Remember to check chromosome before simply incrementing coordinates
	if (wi->chrom == iter->chrom && iter->start < wi->finish + width)
		wi->start = wi->finish;
	else
		// This odd looking formula is an integer division that rounds up
		wi->start = ((iter->start - 1) / width) * width + 1;
	wi->chrom = iter->chrom;
	wi->finish = wi->start + width;

	// Compute sum
	wi->value = 0;
	int total_covered_length = 0;
	while (!iter->done && iter->chrom == wi->chrom && iter->start < wi->finish) {
		int start, finish;
		if (iter->start < wi->start)
			start = wi->start;
		else
			start = iter->start;

		if (iter->finish > wi->finish)
			finish = wi->finish;
		else
			finish = iter->finish;

		int covered_length = finish - start;
		wi->value += covered_length * iter->value;
		total_covered_length += covered_length;
		if (iter->finish > wi->finish)
			break;
		else
			pop(iter);
	}

	// Adding default values if required
	if (iter->default_value && total_covered_length < width)
		wi->value += (width - total_covered_length) * iter->default_value;
}

WiggleIterator * BinningWiggleIterator(WiggleIterator * i, int width) {
	BinningWiggleIteratorData * data = (BinningWiggleIteratorData *) calloc(1, sizeof(BinningWiggleIteratorData));
	if (width < 2) {
		fprintf(stderr, "Cannot bin over a window of width %i, must be 2 or more\n", width);
		exit(1);
	}
	data->iter = NonOverlappingWiggleIterator(i);
	data->width = width;
	return newWiggleIterator(data, &BinningWiggleIteratorPop, &BinningWiggleIteratorSeek, i->default_value * width, false);
}

//////////////////////////////////////////////////////
// Smooth' operator !
//////////////////////////////////////////////////////

typedef struct SmoothWiggleIteratorData_st {
	WiggleIterator * iter;
	double * buffer;
	double sum;
	int latest;
	int oldest;
	int last_position;
	int count;
	int width;
} SmoothWiggleIteratorData;

static double smoothWiggleIteratorRecomputeSum(SmoothWiggleIteratorData * data) {
       double sum = 0;
       int index;
       if (data->oldest < data->latest)
               for (index = data->oldest + 1; index < data->latest; index++)
                       sum += data->buffer[index];
       else {
               for (index = data->oldest + 1; index < data->width; index++)
                       sum += data->buffer[index];
               for (index = 0; index < data->latest; index++)
                       sum += data->buffer[index];
       }

       return sum;
}

static void smoothWiggleIteratorEraseOne(SmoothWiggleIteratorData * data) {
	if (isnan(data->buffer[data->oldest]))
               data->sum = smoothWiggleIteratorRecomputeSum(data);
	else
               data->sum -= data->buffer[data->oldest];
	data->count--;
	if (++data->oldest >= data->width)
		data->oldest = 0;
}

static void smoothWiggleIteratorReadOne(char * chrom, int position, SmoothWiggleIteratorData * data) {
	WiggleIterator * iter = data->iter;

	// Let iter run if necessary
	while (!iter->done && iter->chrom == chrom && iter->finish <= position)
		pop(iter);

	// Record iter's value as appropriate
	if (!iter->done && iter->chrom == chrom) {
		if (iter->start <= position) {
			data->buffer[data->latest] = iter->value;
			data->sum += iter->value;
		} else
			data->buffer[data->latest] = 0;
		data->count++;
		data->last_position = position;

		// Increment ptr
		if (++data->latest >= data->width)
			data->latest = 0;
	}
}

static void SmoothWiggleIteratorPop(WiggleIterator * wi) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done && data->count == 0) {
		// Source is empty, buffer ran out, going home
		wi->done = true;
		return;
	}

	if (data->count == 0) {
		// Smoothing window came to an end
		// Jump to new location
		wi->chrom = iter->chrom;
		wi->start = iter->start - data->width/2;

		if (wi->start < 1) {
			// If start is too close to the start, then we start at zero and prefill the buffer
			wi->start = 1;
			int i;
			for (i = 0; i < data->width/2; i++)
				smoothWiggleIteratorReadOne(wi->chrom, wi->start + i, data);
		}
	} else
		wi->start++;

	// Step by one
	wi->finish = wi->start + 1;
	smoothWiggleIteratorReadOne(wi->chrom, wi->start + data->width/2, data);
	wi->value = data->sum / data->width;

	// Discard unecessary value from buffer
	if (data->count == data->width || data->last_position < wi->start + data->width/2)
		smoothWiggleIteratorEraseOne(data);
}

void SmoothWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	int i;
	for (i = 0; i < data->width; i++)
		data->buffer[i] = 0;
	data->latest = 0;
	data->oldest = 0;
	data->sum = 0;
	data->count = 0;
	wi->done = false;
	pop(wi);
}

WiggleIterator * SmoothWiggleIterator(WiggleIterator * i, int width) {
	SmoothWiggleIteratorData * data = (SmoothWiggleIteratorData *) calloc(1, sizeof(SmoothWiggleIteratorData));
	if (width < 2) {
		fprintf(stderr, "Cannot smooth over a window of width %i, must be 2 or more\n", width);
		exit(1);
	}
	data->iter = NonOverlappingWiggleIterator(i);
	data->buffer = (double*) calloc(sizeof(double), width);
	data->width = width;
	return newWiggleIterator(data, &SmoothWiggleIteratorPop, &SmoothWiggleIteratorSeek, i->default_value, false);
}

//////////////////////////////////////////////////////
// Convenience file reader
//////////////////////////////////////////////////////

WiggleIterator * SmartReader(char * filename, bool holdFire) {
	size_t length = strlen(filename);
	if (!strcmp(filename + length - 3, ".bw"))
		return BigWiggleReader(filename, holdFire);
	else if (!strcmp(filename + length - 7, ".bigWig"))
		return BigWiggleReader(filename, holdFire);
	else if (!strcmp(filename + length - 7, ".bigwig"))
		return BigWiggleReader(filename, holdFire);
	else if (!strcmp(filename + length - 3, ".bg"))
		return WiggleReader(filename);
	else if (!strcmp(filename + length - 4, ".wig"))
		return WiggleReader(filename);
	else if (!strcmp(filename + length - 4, ".bed"))
		return BedReader(filename);
	else if (!strcmp(filename + length - 3, ".bb"))
		return BigBedReader(filename, holdFire);
	else if (!strcmp(filename + length - 7, ".bigBed"))
		return BigBedReader(filename, holdFire);
	else if (!strcmp(filename + length - 7, ".bigbed"))
		return BigBedReader(filename, holdFire);
	else if (!strcmp(filename + length - 4, ".bam"))
		return BamReader(filename, holdFire, false);
	else if (!strcmp(filename + length - 5, ".cram"))
		return BamReader(filename, holdFire, false);
	else if (!strcmp(filename + length - 4, ".sam"))
		return SamReader(filename, false);
	else if (!strcmp(filename + length - 4, ".vcf"))
		return VcfReader(filename);
	else if (!strcmp(filename + length - 4, ".bcf"))
		return BcfReader(filename, holdFire);
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
			iter = data->iter = SmartReader(data->filenames[data->index], false);
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
	data->iter = SmartReader(data->filenames[0], false);
	return newWiggleIterator(data, &CatWiggleIteratorPop, &CatWiggleIteratorSeek, 0, true);
}
