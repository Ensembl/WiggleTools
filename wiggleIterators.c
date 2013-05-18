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

WiggleIterator * newWiggleIterator(void * data, void (*pop)(WiggleIterator *)) {
	WiggleIterator * new = (WiggleIterator *) calloc(1, sizeof(WiggleIterator));
	new->data = data;
	new->pop = pop;
	new->chrom = calloc(1000,1);
	pop(new);
	return new;
}

void destroyWiggleIterator(WiggleIterator * wi) {
	free(wi->data);
	free(wi->chrom);
	free(wi);
}

void pop(WiggleIterator * wi) {
	(*(wi->pop))(wi);
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
		pop(data->iter);
	} else {
		wi->done = true;
	}
}

WiggleIterator * UnitWiggleIterator(WiggleIterator * i) {
	UnitWiggleIteratorData * data = (UnitWiggleIteratorData *) calloc(1, sizeof(UnitWiggleIteratorData));
	data->iter = i;
	return newWiggleIterator(data, &UnitWiggleIteratorPop);
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

WiggleIterator * ScaleWiggleIterator(WiggleIterator * i, double s) {
	ScaleWiggleIteratorData * data = (ScaleWiggleIteratorData *) calloc(1, sizeof(ScaleWiggleIteratorData));
	data->iter = i;
	data->scalar = s;
	return newWiggleIterator(data, &ScaleWiggleIteratorPop);
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

WiggleIterator * NaturalLogWiggleIterator(WiggleIterator * i) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) calloc(1, sizeof(LogWiggleIteratorData));
	data->iter = i;
	data->base = E;
	data->baseLog = 1;
	return newWiggleIterator(data, &LogWiggleIteratorPop);
}

WiggleIterator * LogWiggleIterator(WiggleIterator * i, double s) {
	LogWiggleIteratorData * data = (LogWiggleIteratorData *) calloc(1, sizeof(LogWiggleIteratorData));
	data->iter = i;
	data->base = s;
	data->baseLog = log(s);
	return newWiggleIterator(data, &LogWiggleIteratorPop);
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

WiggleIterator * ExpWiggleIterator(WiggleIterator * i, double s) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) calloc(1, sizeof(ExpWiggleIteratorData));
	data->iter = i;
	data->radix = s;
	data->radixLog = log(data->radix);
	return newWiggleIterator(data, &ExpWiggleIteratorPop);
}

WiggleIterator * NaturalExpWiggleIterator(WiggleIterator * i) {
	ExpWiggleIteratorData * data = (ExpWiggleIteratorData *) calloc(1, sizeof(ExpWiggleIteratorData));
	data->iter = i;
	data->radix = E;
	data->radixLog = 1;
	return newWiggleIterator(data, &ExpWiggleIteratorPop);
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
	return newWiggleIterator(data, &PowerWiggleIteratorPop);
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
	return newWiggleIterator(data, &AbsWiggleIteratorPop);
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
		wi->chrom = iter->chrom;
		wi->finish = -1;
	}
	wi->start = iter->start > wi->finish? iter->start: wi->finish;
	wi->finish = iter->finish;
	wi->value = iter->value;
	pop(iter);
}

void SumWiggleIteratorPop(WiggleIterator * wi) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) wi->data;
	WiggleIterator * iterA = data->iterA;
	WiggleIterator * iterB = data->iterB;
	if (iterA->done && iterB->done) 
		// All done
		wi->done = true;
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
			// Both iterators on the same wi->chromosome:	
			if (!wi->chrom || strcmp(wi->chrom, iterA->chrom) < 0) {
				wi->chrom = iterA->chrom;
				wi->finish = -1;
			}
			if (iterA->start < iterB->start) {
				wi->start = iterA->start > wi->finish? iterA->start: wi->finish;	
				if (iterA->finish <= iterB->start) {
					wi->finish = iterA->finish;
					wi->value = iterA->value;
					pop(iterA);
				} else {
					wi->finish = iterB->start - 1;
					wi->value = iterA->value;
				}
			} else if (iterB->start < iterA->start) {
				wi->start = iterB->start > wi->finish? iterB->start: wi->finish;	
				if (iterB->finish <= iterA->start) {
					wi->finish = iterB->finish;
					wi->value = iterB->value;
					pop(iterB);
				} else {
					wi->finish = iterA->start - 1;
					wi->value = iterB->value;
				}
			} else {
				wi->start = iterA->start > wi->finish? iterA->start: wi->finish;	
				wi->value = iterA->value + iterB->value;
				if (iterA->finish < iterB->finish) {
					pop(iterA);
				} else if (iterB->finish < iterA->finish) {
					pop(iterB);
				} else {
					pop(iterA);
					pop(iterB);
				}
			}
		}
	}
}

WiggleIterator * SumWiggleIterator(WiggleIterator * a, WiggleIterator * b) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) calloc(1, sizeof(BinaryWiggleIteratorData));
	data->iterA = a;
	data->iterB = b;
	return newWiggleIterator(data, &SumWiggleIteratorPop);
}

//////////////////////////////////////////////////////
// Binary Product operator
//////////////////////////////////////////////////////

void ProductWiggleIteratorPop(WiggleIterator * wi) {
	BinaryWiggleIteratorData * data = (BinaryWiggleIteratorData *) wi->data;
	WiggleIterator * iterA = data->iterA;
	WiggleIterator * iterB = data->iterB;

	if (iterA->done || iterB->done)
		wi->done = true;
	else {	
		int chromDiff = strcmp(iterA->chrom, iterB->chrom);
		if (chromDiff < 0) 
			pop(iterA);
		else if (chromDiff > 0)
			pop(iterB);
		else {
			// Both iterators on the same chromosome:	
			if (!wi->chrom || strcmp(wi->chrom, iterA->chrom) < 0) {
				wi->chrom = iterA->chrom;
				wi->finish = -1;
			}
			if (iterA->start < iterB->start) {
				wi->start = iterA->start > wi->finish? iterA->start: wi->finish;	
				if (iterA->finish <= iterB->start) {
					wi->finish = iterA->finish;
					wi->value = 0;
					pop(iterA);
				} else {
					wi->finish = iterB->start - 1;
					wi->value = 0;
				}
			} else if (iterB->start < iterA->start) {
				wi->start = iterB->start > wi->finish? iterB->start: wi->finish;	
				if (iterB->finish <= iterA->start) {
					wi->finish = iterB->finish;
					wi->value = 0;
					pop(iterB);
				} else {
					wi->finish = iterA->start - 1;
					wi->value = 0;
				}
			} else {
				wi->start = iterA->start > wi->finish? iterA->start: wi->finish;	
				wi->value = iterA->value * iterB->value;
				if (iterA->finish < iterB->finish) {
					pop(iterA);
				} else if (iterB->finish < iterA->finish) {
					pop(iterB);
				} else {
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
	return newWiggleIterator(data, &ProductWiggleIteratorPop);
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

WiggleIterator * WigOrBigWigReader(char * filename) {
	size_t length = strlen(filename);
	if (!strcmp(filename + length - 3, ".bw"))
		return BigWiggleReader(filename);
	else if (!strcmp(filename + length - 3, ".bg"))
		return WiggleReader(filename);
	else if (!strcmp(filename + length - 4, ".wig"))
		return WiggleReader(filename);
	else if (!strcmp(filename, "-"))
		return WiggleReader(filename);
	else {
		printf("Could not recognize file format from suffix: %s\n", filename);
		exit(1);
	}
}
