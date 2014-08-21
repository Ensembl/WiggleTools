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
//

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "multiplexer.h"

void popMultiplexer(Multiplexer * multi) {
	if (!multi->done)
		multi->pop(multi);
}

void seekMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	if (!multi->done)
		multi->seek(multi, chrom, start, finish);
}

static void chooseCoords(Multiplexer * multi) {
	int i;
	char * lastChrom = multi->chrom;
	int lastFinish = multi->finish;
	int clipping = -1;
	int comparison;
	int first = -1;
	bool * inplayPtr = multi->inplay;
	WiggleIterator ** wiPtr = multi->iters;

	for (i = 0; i < multi->count; i++) {
		if ((*wiPtr)->done) {
			*inplayPtr = false;
		} else if (first < 0 || (comparison = strcmp((*wiPtr)->chrom, multi->chrom)) < 0) {
			multi->chrom = (*wiPtr)->chrom;

			if (lastChrom && strcmp(multi->chrom, lastChrom) == 0)
				clipping = lastFinish;

			if ((*wiPtr)->start < clipping)
				multi->start = clipping;
			else
				multi->start = (*wiPtr)->start;

			multi->finish = (*wiPtr)->finish;
			*inplayPtr = true;
			first = i;
		} else if (comparison > 0){
			*inplayPtr = false;
		} else {
			int start = (*wiPtr)->start;
			if (start > multi->start) {
				*inplayPtr = false;
				if (start < multi->finish)
					multi->finish = start;
			} else {

				*inplayPtr = true;

				if (multi->start > clipping && start < multi->start) {
					multi->finish = multi->start;
					if (start < clipping)
						multi->start = clipping;
					else
						multi->start = start;
					first = i;
				}

				if ((*wiPtr)->finish < multi->finish)
					multi->finish = (*wiPtr)->finish;
			}
		}

		inplayPtr++;
		wiPtr++;
	}

	if (first < 0) {
		multi->done = true;
		return;
	}

	for (i = 0; i < first; i++)
		multi->inplay[i] = false;
}

static bool readValues(Multiplexer * multi) {
	int i;
	bool * inplayPtr = multi->inplay;
	WiggleIterator ** wiPtr = multi->iters;
	double * valuePtr = multi->values;
	bool allActive = true;

	for (i=0; i < multi->count; i++) {
		if (multi->strict && (*wiPtr)->done) {
			multi->done;
			break;
		}

		if (*inplayPtr) {
			*valuePtr = (*wiPtr)->value;
			if ((*wiPtr)->finish == multi->finish) {
				pop(*wiPtr);
			}
		} else 
			allActive = false;

		valuePtr++;
		inplayPtr++;
		wiPtr++;
	}

	return allActive;
}

static void popCoreMultiplexer(Multiplexer * multi) {
	while (!multi->done) {
		chooseCoords(multi);
		if (multi->done)
			break;
		if (readValues(multi) || !multi->strict)
			break;
	}
}

static void seekCoreMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	int i;
	multi->done = false;
	for (i=0; i<multi->count; i++)
		seek(multi->iters[i], chrom, start, finish);
	multi->chrom = NULL;
	multi->start = 0;
	popMultiplexer(multi);
}

Multiplexer * newCoreMultiplexer(WiggleIterator ** iters, int count, void (*pop)(Multiplexer *), void (*seek)(Multiplexer *, const char *, int, int)) {
	Multiplexer * new = (Multiplexer *) calloc (1, sizeof(Multiplexer));
	new->count = count;
	new->iters = calloc(count, sizeof(WiggleIterator *));
	int i;
	for (i = 0; i < count; i++)
		new->iters[i] = NonOverlappingWiggleIterator(iters[i]);
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double *) calloc(count, sizeof(double));
	new->pop = pop;
	new->seek = seek;
	return new;
}

Multiplexer * newMultiplexer(WiggleIterator ** iters, int count, bool strict) {
	Multiplexer * new = newCoreMultiplexer(iters, count, popCoreMultiplexer, seekCoreMultiplexer);
	new->strict = strict;
	popMultiplexer(new);
	return new;
}
