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

#include "multiSet.h"

static void popClosingMultiplexers(Multiset * multi) {
	while (fh_notempty(multi->finishes) && fh_min(multi->finishes) == multi->finish) {
		int index = fh_extractmin(multi->finishes);
		Multiplexer * multiplexer = multi->multis[index];
		popMultiplexer(multiplexer);
		multi->inplay[index] = false;
		multi->inplay_count--;
		if (!multiplexer->done && !strcmp(multiplexer->chrom, multi->chrom))
			fh_insert(multi->starts, multiplexer->start, index);
	}
}

static void queueUpMultiplexers(Multiset * multi) {
	// Find lowest value chromosome
	multi->chrom = NULL;
	Multiplexer ** muPtr = multi->multis;
	int i; 
	for (i = 0; i < multi->count; i++) {
		if ((!(*muPtr)->done) && (!multi->chrom || strcmp((*muPtr)->chrom, multi->chrom) < 0))
			multi->chrom = (*muPtr)->chrom;
		muPtr++;
	}

	// No chromosome found => All multisets done
	if (!multi->chrom) {
		multi->done = true;
		return;
	}

	// Put those multiplexers in heap
	muPtr = multi->multis;
	for (i = 0; i < multi->count; i++) {
		if ((!(*muPtr)->done) && strcmp((*muPtr)->chrom, multi->chrom) == 0)
			fh_insert(multi->starts, (*muPtr)->start, i);
		muPtr++;
	}

}

static void admitNewMultiplexersIntoPlay(Multiset * multi) {
	while (fh_notempty(multi->starts) && fh_min(multi->starts) == multi->start) {
		int index = fh_extractmin(multi->starts);
		Multiplexer * multiplexer = multi->multis[index];
		fh_insert(multi->finishes, multiplexer->finish, index);
		multi->inplay[index] = true;
		multi->inplay_count++;
	}
}

static void defineNewFinish(Multiset * multi) {
	multi->finish = fh_min(multi->finishes);

	if (fh_notempty(multi->starts)) {
		int min_start = fh_min(multi->starts);
		if (multi->finish > min_start)
			multi->finish = min_start;
	}
}

void popMultiset(Multiset * multi) {
	popClosingMultiplexers(multi);

	// Check that there are multiplexers queued up
	// If no multiplexers are waiting, either waiting on other chromosomes
	// or finished.
	if (fh_empty(multi->starts) && fh_empty(multi->finishes))
		queueUpMultiplexers(multi);

	// If queues still empty
	if (multi->done)
		return;

	// If no multiplexer in play jump to next start
	if (multi->inplay_count)
		multi->start = multi->finish;
	else
		multi->start = fh_min(multi->starts);

	admitNewMultiplexersIntoPlay(multi);
	defineNewFinish(multi);
}

void seekMultiset(Multiset * multi, const char * chrom, int start, int finish) {
	int i;
	multi->done = false;
	for (i=0; i<multi->count; i++)
		seekMultiplexer(multi->multis[i], chrom, start, finish);
	fh_deleteheap(multi->starts);
	fh_deleteheap(multi->finishes);
	multi->starts = fh_makeheap();
	multi->finishes = fh_makeheap();
	popMultiset(multi);
}

Multiset * newMultiset(Multiplexer ** multis, int count) {
	Multiset * new = (Multiset *) calloc (1, sizeof(Multiset));
	new->count = count;
	new->multis = multis;
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double **) calloc(count, sizeof(double));
	new->starts = fh_makeheap();
	new->finishes = fh_makeheap();
	int i;
	for (i = 0; i < count; i++)
		new->values[i] = multis[i]->values;
	popMultiset(new);
	return new;
}
