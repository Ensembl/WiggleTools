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
//

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "multiplexer.h"

void popMultiplexer(Multiplexer * multi) {
	if (!multi->done)
		multi->pop(multi);
}

void runMultiplexer(Multiplexer * multi) {
	while (!multi->done)
		multi->pop(multi);
}

void seekMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	multi->done = false;
	multi->seek(multi, chrom, start, finish);
}

static void popClosingWiggleIterators(Multiplexer * multi) {
	while (fh_notempty(multi->finishes) && fh_min(multi->finishes) == multi->finish) {
		int index = fh_extractmin(multi->finishes);
		WiggleIterator * wi = multi->iters[index];
		pop(wi);
		multi->inplay[index] = false;
		multi->inplay_count--;
		multi->values[index] = wi->default_value;
		if (!wi->done && !strcmp(wi->chrom, multi->chrom))
			fh_insert(multi->starts, wi->start, index);
	}
}

static void queueUpWiggleIterators(Multiplexer * multi) {
	// Find lowest value chromosome
	multi->chrom = NULL;
	WiggleIterator ** muPtr = multi->iters;
	int i; 
	for (i = 0; i < multi->count; i++) {
		if ((!(*muPtr)->done) && (!multi->chrom || strcmp((*muPtr)->chrom, multi->chrom) < 0))
			multi->chrom = (*muPtr)->chrom;
		muPtr++;
	}

	// No chromosome found => All itersets done
	if (!multi->chrom) {
		multi->done = true;
		return;
	}

	// Put those wis in heap
	muPtr = multi->iters;
	for (i = 0; i < multi->count; i++) {
		if ((!(*muPtr)->done) && strcmp((*muPtr)->chrom, multi->chrom) == 0)
			fh_insert(multi->starts, (*muPtr)->start, i);
		muPtr++;
	}
}

static void admitNewWiggleIteratorsIntoPlay(Multiplexer * multi) {
	while (fh_notempty(multi->starts) && fh_min(multi->starts) == multi->start) {
		int index = fh_extractmin(multi->starts);
		WiggleIterator * wi = multi->iters[index];
		fh_insert(multi->finishes, wi->finish, index);
		multi->inplay[index] = true;
		multi->values[index] = wi->value;
		multi->inplay_count++;
	}
}

static void defineNewFinish(Multiplexer * multi) {
	multi->finish = fh_min(multi->finishes);

	if (fh_notempty(multi->starts)) {
		int min_start = fh_min(multi->starts);
		if (multi->finish > min_start) {
			multi->finish = min_start;
		}
	}
}

static bool popCoreMultiplexer2(Multiplexer * multi) {
	popClosingWiggleIterators(multi);

	// Check that there are wis queued up
	// If no wis are waiting, either waiting on other chromosomes
	// or finished.
	if (fh_empty(multi->starts) && fh_empty(multi->finishes))
		queueUpWiggleIterators(multi);

	// If queues still empty
	if (multi->done)
		return false;

	// If no wi in play jump to next start
	if (multi->inplay_count)
		multi->start = multi->finish;
	else
		multi->start = fh_min(multi->starts);

	admitNewWiggleIteratorsIntoPlay(multi);
	defineNewFinish(multi);

	return multi->inplay_count == multi->count;
}

static void popCoreMultiplexer(Multiplexer * multi) {
	while (!multi->done) {
		if (popCoreMultiplexer2(multi) || !multi->strict)
			break;
	}
}

static void seekCoreMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	int i;
	multi->done = false;
	for (i=0; i<multi->count; i++)
		seek(multi->iters[i], chrom, start, finish);
	fh_deleteheap(multi->starts);
	fh_deleteheap(multi->finishes);
	multi->starts = fh_makeheap();
	multi->finishes = fh_makeheap();
	multi->inplay_count = 0;
	popMultiplexer(multi);
}

Multiplexer * newCoreMultiplexer(void * data, int count, void (*pop)(Multiplexer *), void (*seek)(Multiplexer *, const char *, int, int)) {
	Multiplexer * new = (Multiplexer *) calloc (1, sizeof(Multiplexer));
	new->count = count;
	new->values = (double *) calloc(count, sizeof(double));
	new->default_values = (double *) calloc(count, sizeof(double));
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->pop = pop;
	new->seek = seek;
	new->data = data;
	new->starts = fh_makeheap();
	new->finishes = fh_makeheap();
	return new;
}

Multiplexer * newMultiplexer(WiggleIterator ** iters, int count, bool strict) {
	Multiplexer * new = newCoreMultiplexer(NULL, count, popCoreMultiplexer, seekCoreMultiplexer);
	new->strict = strict;
	new->iters = calloc(count, sizeof(WiggleIterator *));
	int i;
	for (i = 0; i < count; i++) {
		new->iters[i] = NonOverlappingWiggleIterator(iters[i]);
		new->default_values[i] = new->iters[i]->default_value;
		new->values[i] = new->iters[i]->default_value;
	}
	popMultiplexer(new);
	return new;
}
