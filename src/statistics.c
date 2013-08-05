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
#include "wiggleIterator.h"

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

