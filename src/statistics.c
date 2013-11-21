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

double isZero(WiggleIterator * wi) {
	for (; !wi->done; pop(wi))
		if (wi->value != 0)
			exit(1);
	return 1;
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
// Profile summaries
//////////////////////////////////////////////////////

void regionProfile2(WiggleIterator * wig, int width, int centerRes, int centerRegion, double scaleFactor, double * res, bool stranded, int strand) {
	int start, finish, pos;

	if (!stranded || strand > 0) {
		start = centerRes + (wig->start - centerRegion); 
		finish = centerRes + (wig->finish - centerRegion); 
	} else if (strand < 0) {
		start = centerRes - (wig->finish - centerRegion); 
		finish = centerRes - (wig->start - centerRegion); 
	} else {
		fprintf(stderr, "Cannot provide stranded profile on non-stranded regions\n");
		exit(1);
	}

	for (pos = start; pos < finish; pos++) {
		int index = (int) round(pos * scaleFactor);
		if (index < 0)
			index = 0;
		if (index >= width)
			index = width - 1;
		res[index] += wig->value;
	}
}

void regionProfile(WiggleIterator * region, WiggleIterator * wig, int width, double * res, int centerRes, bool stranded) {
	int centerRegion = (region->finish + region->start) / 2;
	double scaleFactor = width / (double) (region->finish - region->start);

	for (seek(wig, region->chrom, region->start, region->finish); !wig->done; pop(wig)) 
		regionProfile2(wig, width, centerRes, centerRegion, scaleFactor, res, stranded, region->strand);
}

double * profileSum(WiggleIterator * regions, WiggleIterator * wig, int width, bool stranded) {
	double * res = calloc(width, sizeof(double));
	int centerRes = width / 2; 

	for (; !regions->done; pop(regions)) 
		regionProfile(regions, wig, width, res, centerRes, stranded);

	return res;
}
