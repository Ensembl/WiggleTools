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
	WiggleIterator * wi =  newWiggleIterator(data, &CompressionWiggleIteratorPop, &UnitWiggleIteratorSeek);
	return wi;
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
// Correlation function
//////////////////////////////////////////////////////
// Note: this is an approximate calculation of the Pearson calculation
// which has the benefit of running in a single pass through the data
// The origin of the code can be found at:
// http://en.wikipedia.org/wiki?title=Talk:Correlation

double pearsonCorrelation(WiggleIterator * iterA, WiggleIterator * iterB) {
	double sum_sq_A = 0;
	double sum_sq_B = 0;
	double sum_AB = 0;
	double meanA = iterA->value;
	double meanB = iterB->value;
	int totalLength = 0;
	int halfway, width;
	double sweep, deltaA, deltaB;
	char *chrom = NULL; 
	int start = -1;
	int finish = -1;

	while (!iterA->done || !iterB->done) {
		if (iterA->done) {
			if (!chrom || strcmp(chrom, iterB->chrom) < 0) {
				chrom = iterB->chrom;
				finish = -1;
			}
			start = iterB->start > finish? iterB->start: finish;
			finish = iterB->finish;
			width = (finish - start);
			halfway = totalLength + width/2;
			if (halfway == 0)
				halfway = 1;
			sweep = (halfway - 1.0) / halfway;
			totalLength += width;
			deltaA = -meanA;
			deltaB = iterB->value - meanB;
			sum_sq_A += deltaA * deltaA * sweep;
			sum_sq_B += deltaB * deltaB * sweep;
			sum_AB += deltaA * deltaB * sweep;
			meanA += deltaA / totalLength;
			meanB += deltaB / totalLength;
			pop(iterB);
			continue;
		}
		if (iterB->done) {

			if (!chrom || strcmp(chrom, iterA->chrom) < 0) {
				chrom = iterA->chrom;
				finish = -1;
			}
			start = iterA->start > finish? iterA->start: finish;
			finish = iterA->finish;
			width = (finish - start);
			halfway = totalLength + width/2;
			if (halfway == 0)
				halfway = 1;
			sweep = (halfway - 1.0) / halfway;
			totalLength += width;
			deltaB = -meanB;
			deltaA = iterA->value - meanA;
			sum_sq_B += deltaB * deltaB * sweep;
			sum_sq_A += deltaA * deltaA * sweep;
			sum_AB += deltaB * deltaA * sweep;
			meanB += deltaB / totalLength;
			meanA += deltaA / totalLength;
			pop(iterA);
			continue;
		}

		int chromDiff = strcmp(iterA->chrom, iterB->chrom);

		if (chromDiff < 0) {
			if (!chrom || strcmp(chrom, iterA->chrom) < 0) {
				chrom = iterA->chrom;
				finish = -1;
			}
			start = iterA->start > finish? iterA->start: finish;
			finish = iterA->finish;
			width = (finish - start);
			halfway = totalLength + width/2;
			if (halfway == 0)
				halfway = 1;
			sweep = (halfway - 1.0) / halfway;
			totalLength += width;
			deltaB = -meanB;
			deltaA = iterA->value - meanA;
			sum_sq_B += deltaB * deltaB * sweep;
			sum_sq_A += deltaA * deltaA * sweep;
			sum_AB += deltaB * deltaA * sweep;
			meanB += deltaB / totalLength;
			meanA += deltaA / totalLength;
			pop(iterA);
		} else if (chromDiff > 0) {
			if (!chrom || strcmp(chrom, iterB->chrom) < 0) {
				chrom = iterB->chrom;
				finish = -1;
			}
			start = iterB->start > finish? iterB->start: finish;
			finish = iterB->finish;
			width = (finish - start);
			halfway = totalLength + width/2;
			if (halfway == 0)
				halfway = 1;
			sweep = (halfway - 1.0) / halfway;
			totalLength += width;
			deltaA = -meanA;
			deltaB = iterB->value - meanB;
			sum_sq_A += deltaA * deltaA * sweep;
			sum_sq_B += deltaB * deltaB * sweep;
			sum_AB += deltaA * deltaB * sweep;
			meanA += deltaA / totalLength;
			meanB += deltaB / totalLength;
			pop(iterB);
		} else {
			// Both iterators on the same chromosome:	
			chrom = iterA->chrom;

			if (iterA->start < iterB->start) {
				start = iterA->start > finish? iterA->start: finish;	
				if (iterA->finish <= iterB->start) {
					finish = iterA->finish;
				} else {
					finish = iterB->start - 1;
				}
				deltaB = -meanB;
				deltaA = iterA->value - meanA;
			} else if (iterB->start < iterA->start) {
				start = iterB->start > finish? iterB->start: finish;	
				if (iterB->finish <= iterA->start) {
					finish = iterB->finish;
				} else {
					finish = iterA->start - 1;
				}
				deltaA = -meanA;
				deltaB = iterB->value - meanB;
			} else {
				start = iterA->start > finish? iterA->start: finish;	
				finish = iterA->finish < iterB->finish ? iterA->finish: iterB->finish;
				deltaA = iterA->value - meanA;
				deltaB = iterB->value - meanB;
			}

			width = (finish - start);
			halfway = totalLength + width/2;
			if (halfway == 0)
				halfway = 1;
			sweep = (halfway - 1.0) / halfway;
			totalLength += width;
			deltaA = iterA->value - meanA;
			deltaB = iterB->value - meanB;
			sum_sq_A += deltaA * deltaA * sweep;
			sum_sq_B += deltaB * deltaB * sweep;
			sum_AB += deltaA * deltaB * sweep;
			meanA += deltaA / totalLength;
			meanB += deltaB / totalLength;

			if (iterA->finish == finish)
				pop(iterA);
			if (iterB->finish == finish)
				pop(iterB);
		}
	}

	return sum_AB / (sqrt(sum_sq_A) * sqrt(sum_sq_B));
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
