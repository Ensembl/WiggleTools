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

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Local header
#include "wiggleTools.h"
#include "wiggleIterators.h"

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
	WiggleIterator * new = newWiggleIterator(NULL, &NullWiggleIteratorPop, &NullWiggleIteratorSeek);
	new->done = true;
	new->done = true;
	return new;
}

//////////////////////////////////////////////////////
// Unit operator
//////////////////////////////////////////////////////

typedef struct UnitWiggleIteratorData_st {
	WiggleIterator * iter;
} UnitWiggleIteratorData;

void UnitWiggleIteratorPop(WiggleIterator * wi) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		if (iter->value > 0)
			wi->value = 1;
		else if (iter->value < 0)
			wi->value = -1;
		else
			wi->value = 0;
		pop(iter);
	} else {
		wi->done = true;
	}
}

void UnitWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
	wi->done = false;
	pop(wi);
}

WiggleIterator * UnitWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &UnitWiggleIteratorPop, &UnitWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Compression operator
//////////////////////////////////////////////////////

void CompressionWiggleIteratorPop(WiggleIterator * wi) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;

	if (iter->done) {
		wi->done = true;
	} else {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;
		pop(iter);

		while (!iter->done && strcmp(iter->chrom, wi->chrom) == 0 && iter->start == wi->finish && iter->value == wi->value) {
			wi->finish = iter->finish;
			pop(iter);
		}
	}
}


WiggleIterator * CompressionWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &CompressionWiggleIteratorPop, &UnitWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Union operator
//////////////////////////////////////////////////////

void UnionWiggleIteratorPop(WiggleIterator * wi) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
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
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	WiggleIterator * wi =  newWiggleIterator(data, &UnionWiggleIteratorPop, &UnitWiggleIteratorSeek);
	wi->value = 1;
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
	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
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
	data->iter = i;
	data->scalar = s;
	return newWiggleIterator(data, &ScaleWiggleIteratorPop, &ScaleWiggleIteratorSeek);
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
	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
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
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek);
}

WiggleIterator * LogWiggleIterator(WiggleIterator * i, double s) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) calloc(1, sizeof(LogWiggleIteratorData));
	data->iter = i;
	data->base = s;
	data->baseLog = log(s);
	return newWiggleIterator(data, &LogWiggleIteratorPop, &LogWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Exponentiation operator
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
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek);
}

WiggleIterator * NaturalExpWiggleIterator(WiggleIterator * i) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) calloc(1, sizeof(ExpWiggleIteratorData));
	data->iter = i;
	data->radix = E;
	data->radixLog = 1;
	return newWiggleIterator(data, &ExpWiggleIteratorPop, &ExpWiggleIteratorSeek);
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
		wi->value = pow(iter->value, data->scalar);
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * PowerWiggleIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = i;
	data->scalar = s;
	return newWiggleIterator(data, &PowerWiggleIteratorPop, &ScaleWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Abs operator
//////////////////////////////////////////////////////

static void AbsWiggleIteratorPop(WiggleIterator * wi) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = abs(iter->value);
		pop(iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * AbsWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &AbsWiggleIteratorPop, &UnitWiggleIteratorSeek);
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
	data->iter = i;
	data->buffer = (double*) calloc(1, width);
	data->width = width;
	return newWiggleIterator(data, &SmoothWiggleIteratorPop, &SmoothWiggleIteratorSeek);
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
	else if (!strcmp(filename + length - 4, ".bin"))
		return BinaryFileReader(filename);
	else if (!strcmp(filename, "-"))
		return WiggleReader(filename);
	else {
		printf("Could not recognize file format from suffix: %s\n", filename);
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
	return newWiggleIterator(data, &CatWiggleIteratorPop, &CatWiggleIteratorSeek);
}
