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
#include <stdio.h>
#include <string.h>

// Local header
#include "wiggleTools.h"
#include "wiggleIterators.h"

void goToNextNonZeroValue(WiggleIterator *  wi) {
	while (!wi->nextDone) {
		(*(wi->pop))(wi);
		if (wi->nextValue)
			break;
	}
}

void oneStep(WiggleIterator * wi) {
	wi->chrom = wi->nextChrom;	
	wi->start = wi->nextStart;
	wi->finish = wi->nextFinish;
	wi->value = wi->nextValue;
	wi->done = wi->nextDone;
}

void pop(WiggleIterator * wi) {
	// Safety check
	if (wi->done)
		return;

	// Update
	oneStep(wi);
	goToNextNonZeroValue(wi);

	// Compress any consecutive blocks
	while (!wi->nextDone && (!strcmp(wi->chrom, wi->nextChrom) && wi->finish == wi->nextStart && wi->value == wi->nextValue)) {
		wi->finish = wi->nextFinish;
		goToNextNonZeroValue(wi);
	}
}

WiggleIterator * newWiggleIterator(void * data, void (*popFunction)(WiggleIterator *), void (*seek)(WiggleIterator *, const char *, int, int)) {
	WiggleIterator * new = (WiggleIterator *) calloc(1, sizeof(WiggleIterator));
	new->data = data;
	new->pop = popFunction;
	new->seek = seek;
	new->chrom = calloc(1000,1);
	new->nextChrom = calloc(1000,1);
	new->nextValue = 1; // Default value for non-valued bed tracks;
	pop(new);
	pop(new);
	return new;
}

void destroyWiggleIterator(WiggleIterator * wi) {
	free(wi->data);
	free(wi->chrom);
	free(wi);
}

void seek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	(*(wi->seek))(wi, chrom, start, finish);
	pop(wi);
}

FILE * openOrFail(char * filename, char * description, char * mode) {
	FILE * file;
	if (!(file = fopen(filename, mode))) {
		printf("Could not open %s %s, exiting...\n", (char *) description, (char *) filename);
	}
	return file;
}

//////////////////////////////////////////////////////
// Output
//////////////////////////////////////////////////////

static void print(WiggleIterator * wi, FILE * out) {
	while (!wi->done) {
		fprintf(out, "%s\t%i\t%i\t%f\n", wi->chrom, wi->start, wi->finish - 1, wi->value);
		pop(wi);
	}
}

void toFile(WiggleIterator * wi, char * filename) {
	FILE * file = openOrFail(filename, "output file", "w");
	print(wi, file);
	fclose(file);
}

void toStdout(WiggleIterator * wi) {
	print(wi, stdout);
}

//////////////////////////////////////////////////////
// Basic Stats
//////////////////////////////////////////////////////

double AUC(WiggleIterator * wi) {
	double total = 0;
	for(;!wi->done; pop(wi)) 
		total += (wi->finish - wi->start) * wi->value;
	return total;
}

double span(WiggleIterator * wi) {
	double total = 0;
	for(;!wi->done; pop(wi)) 
		if (wi->value)
			total += (wi->finish - wi->start);
	return total;
}

double mean(WiggleIterator * wi) {
	double total = 0;
	double span = 0;
	for(;!wi->done; pop(wi)) {
		span += (wi->finish - wi->start);
		total += (wi->finish - wi->start) * wi->value;
	}
	return total / span;
}

double variance(WiggleIterator * wi) {
	// Online algorithm copied from 
	// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass_algorithm
	double sumWeight = 0;
	double mean = 0;
	double M2 = 0;
	double count = 0;

	for(;!wi->done; pop(wi)) {
		double weight = wi->finish - wi->start;
		double temp = sumWeight + weight;
		double delta = wi->value - mean;
		double R = delta * weight / temp;
		mean += R;
		M2 += sumWeight * delta * R;
		sumWeight = temp;
		count++;
	}
	return (M2 * count) / (sumWeight * (count - 1));
}

double stddev(WiggleIterator * wi) {
	return sqrt(variance(wi));
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		if (iter->value > 0)
			wi->nextValue = 1;
		else if (iter->value < 0)
			wi->nextValue = -1;
		else
			wi->nextValue = 0;
		pop(iter);
	} else {
		wi->nextDone = true;
	}
}

void UnitWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
}

WiggleIterator * UnitWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &UnitWiggleIteratorPop, &UnitWiggleIteratorSeek);
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		wi->nextValue = data->scalar * iter->value;
		pop(data->iter);
	} else {
		wi->nextDone = true;
	}
}

void ScaleWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		wi->nextValue = log(iter->value) / data->baseLog;
		pop(data->iter);
	} else {
		wi->nextDone = true;
	}
}

void LogWiggleIteratorSeek(WiggleIterator * wi, const char  * chrom, int start, int finish) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		wi->nextValue = exp(iter->value * data->radixLog);
		pop(iter);
	} else {
		wi->nextDone = true;
	}
}

void ExpWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) wi->data;
	seek(data->iter, chrom, start, finish);
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		wi->nextValue = pow(iter->value, data->scalar);
		pop(iter);
	} else {
		wi->nextDone = true;
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
		wi->nextChrom = iter->chrom;
		wi->nextStart = iter->start;
		wi->nextFinish = iter->finish;
		wi->nextValue = abs(iter->value);
		pop(iter);
	} else {
		wi->nextDone = true;
	}
}

WiggleIterator * AbsWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &AbsWiggleIteratorPop, &UnitWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Binary Sum operator
//////////////////////////////////////////////////////

typedef struct BinaryWiggleIteratorData_st {
	WiggleIterator * iterA;
	WiggleIterator * iterB;
} BinaryWiggleIteratorData;

static void copyInterval(WiggleIterator * wi, WiggleIterator * iter) {
	if (!wi->chrom || strcmp(wi->chrom, iter->chrom) < 0) {
		wi->nextChrom = iter->chrom;
		wi->nextFinish = -1;
	}
	wi->nextStart = iter->start > wi->nextFinish? iter->start: wi->nextFinish;
	wi->nextFinish = iter->finish;
	wi->nextValue = iter->value;
	pop(iter);
}

void SumWiggleIteratorPop(WiggleIterator * wi) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) wi->data;
	WiggleIterator * iterA = data->iterA;
	WiggleIterator * iterB = data->iterB;
	if (iterA->done && iterB->done) 
		// All done
		wi->nextDone = true;
	else if (iterA->done)
		// A expired
		copyInterval(wi, iterB);
	else if (iterB->done)
		// B expired
		copyInterval(wi, iterA);
	else {
		int chromDiff = strcmp(iterA->chrom, iterB->chrom);

		if (chromDiff < 0) 
			// A on previous chromosome
			copyInterval(wi, iterA);
		else if (chromDiff > 0)
			// B on previous chromosome
			copyInterval(wi, iterB);
		else {
			// Both iterators on the same wi->nextChromosome:	
			if (!wi->nextChrom || strcmp(wi->nextChrom, iterA->chrom) < 0) {
				wi->nextChrom = iterA->chrom;
				wi->nextFinish = -1;
			}
			if (iterA->start < iterB->start) {
				wi->nextStart = iterA->start > wi->nextFinish? iterA->start: wi->nextFinish;	
				if (iterA->finish <= iterB->start) {
					wi->nextFinish = iterA->finish;
					wi->nextValue = iterA->value;
					pop(iterA);
				} else {
					wi->nextFinish = iterB->start - 1;
					wi->nextValue = iterA->value;
				}
			} else if (iterB->start < iterA->start) {
				wi->nextStart = iterB->start > wi->nextFinish? iterB->start: wi->nextFinish;	
				if (iterB->finish <= iterA->start) {
					wi->nextFinish = iterB->finish;
					wi->nextValue = iterB->value;
					pop(iterB);
				} else {
					wi->nextFinish = iterA->start - 1;
					wi->nextValue = iterB->value;
				}
			} else {
				wi->nextStart = iterA->start > wi->nextFinish? iterA->start: wi->nextFinish;	
				wi->nextValue = iterA->value + iterB->value;
				if (iterA->finish < iterB->finish) {
					wi->nextFinish = iterA->finish;
					pop(iterA);
				} else if (iterB->finish < iterA->finish) {
					wi->nextFinish = iterB->finish;
					pop(iterB);
				} else {
					wi->nextFinish = iterA->finish;
					pop(iterA);
					pop(iterB);
				}
			}
		}
	}
}

void BinaryWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) wi->data;
	seek(data->iterA, chrom, start, finish);
	seek(data->iterB, chrom, start, finish);

}
WiggleIterator * SumWiggleIterator(WiggleIterator * a, WiggleIterator * b) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) calloc(1, sizeof(BinaryWiggleIteratorData));
	data->iterA = a;
	data->iterB = b;
	return newWiggleIterator(data, &SumWiggleIteratorPop, &BinaryWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Binary Product operator
//////////////////////////////////////////////////////

void ProductWiggleIteratorPop(WiggleIterator * wi) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) wi->data;
	WiggleIterator * iterA = data->iterA;
	WiggleIterator * iterB = data->iterB;

	if (iterA->done || iterB->done)
		wi->nextDone = true;
	else {	
		int chromDiff = strcmp(iterA->chrom, iterB->chrom);
		if (chromDiff < 0) 
			pop(iterA);
		else if (chromDiff > 0)
			pop(iterB);
		else {
			// Both iterators on the same chromosome:	
			if (!wi->nextChrom || strcmp(wi->nextChrom, iterA->chrom) < 0) {
				wi->nextChrom = iterA->chrom;
				wi->nextFinish = -1;
			}
			if (iterA->start < iterB->start) {
				wi->nextStart = iterA->start > wi->nextFinish? iterA->start: wi->nextFinish;	
				if (iterA->finish <= iterB->start) {
					wi->nextFinish = iterA->finish;
					wi->nextValue = 0;
					pop(iterA);
				} else {
					wi->nextFinish = iterB->start - 1;
					wi->nextValue = 0;
				}
			} else if (iterB->start < iterA->start) {
				wi->nextStart = iterB->start > wi->nextFinish? iterB->start: wi->nextFinish;	
				if (iterB->finish <= iterA->start) {
					wi->nextFinish = iterB->finish;
					wi->nextValue = 0;
					pop(iterB);
				} else {
					wi->nextFinish = iterA->start - 1;
					wi->nextValue = 0;
				}
			} else {
				wi->nextStart = iterA->start > wi->nextFinish? iterA->start: wi->nextFinish;	
				wi->nextValue = iterA->value * iterB->value;
				if (iterA->finish < iterB->finish) {
					wi->nextFinish = iterA->finish;
					pop(iterA);
				} else if (iterB->finish < iterA->finish) {
					wi->nextFinish = iterB->finish;
					pop(iterB);
				} else {
					wi->nextFinish = iterA->finish;
					pop(iterA);
					pop(iterB);
				}
			}
		}
	}
}

WiggleIterator * ProductWiggleIterator(WiggleIterator * a, WiggleIterator * b) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) calloc(1, sizeof(BinaryWiggleIteratorData));
	data->iterA = a;
	data->iterB = b;
	return newWiggleIterator(data, &ProductWiggleIteratorPop, &BinaryWiggleIteratorSeek);
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
// Regional stat iterator:
//////////////////////////////////////////////////////

typedef struct applyWiggleIteratorData_st {
	WiggleIterator * regions;
	double (*statistic)(WiggleIterator *);
	WiggleIterator * data;
} ApplyWiggleIteratorData;

void ApplyWiggleIteratorPop(WiggleIterator * wi) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) wi->data;
	wi->chrom = data->regions->chrom;
	wi->start = data->regions->start;
	wi->finish = data->regions->finish;
	seek(wi->data, wi->chrom, wi->start, wi->finish);
	wi->value = (*(data->statistic))(data->data);
}

void ApplyWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) wi->data;
	seek(data->regions, chrom, start, finish);
}

WiggleIterator * apply(WiggleIterator * regions, double (*statistic)(WiggleIterator *), WiggleIterator * dataset) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) calloc(1, sizeof(ApplyWiggleIteratorData));
	data->regions = regions;
	data->statistic = statistic;
	data->data = dataset;
	return newWiggleIterator(data, &ApplyWiggleIteratorPop, &ApplyWiggleIteratorSeek);
}

//////////////////////////////////////////////////////
// Convenience file reader
//////////////////////////////////////////////////////

WiggleIterator * WigOrBigWigReader(char * filename) {
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
	else if (!strcmp(filename, "-"))
		return WiggleReader(filename);
	else {
		printf("Could not recognize file format from suffix: %s\n", filename);
		exit(1);
	}
}
