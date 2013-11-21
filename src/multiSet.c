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

#include "multiSet.h"

void popMultiset(Multiset * multi) {
	(*(multi->pop))(multi);
}

static void chooseCoords(Multiset * multi) {
	int i; 
	char * lastChrom = multi->chrom;
	int lastFinish = multi->finish;
	int clipping = -1;
	int comparison;
	int first = -1;
	bool * inplayPtr = multi->inplay;
	Multiplexer ** muPtr = multi->multis;

	for (i = 0; i < multi->count; i++) {
		if ((*muPtr)->done) {
			*inplayPtr = false;
		} else if (first < 0 || (comparison = strcmp((*muPtr)->chrom, multi->chrom)) < 0) {
			multi->chrom = (*muPtr)->chrom;

			if (lastChrom && strcmp(multi->chrom, lastChrom) == 0)
				clipping = lastFinish;

			if ((*muPtr)->start < clipping)
				multi->start = clipping;
			else
				multi->start = (*muPtr)->start;

			multi->finish = (*muPtr)->finish;
			*inplayPtr = true;
			first = i;
		} else if (comparison > 0){
			*inplayPtr = false;
		} else {
			int start = (*muPtr)->start;
			*inplayPtr = true;
			if (multi->start > clipping && start < multi->start) {
				multi->finish = multi->start;
				if (start < clipping)
					multi->start = clipping;
				else
					multi->start = start;
				first = i;
			} else if (start > multi->start) {
				*inplayPtr = false;
			} 
			
			if ((*muPtr)->finish < multi->finish) {
				multi->finish = (*muPtr)->finish;
			} 
		}

		inplayPtr++;
		muPtr++;
	}

	if (first < 0) {
		multi->done = true;
		return;
	}

	for (i = 0; i < first; i++) 
		multi->inplay[i] = false;

}

void popMultiplexers(Multiset * multi) {
	int i; 
	bool * inplayPtr = multi->inplay;
	Multiplexer ** mpPtr = multi->multis;

	for (i=0; i < multi->count; i++) {
		if (*inplayPtr) {
			if ((*mpPtr)->finish == multi->finish) {
				popMultiplexer(*mpPtr);
			}
		}
		inplayPtr++;
		mpPtr++;
	}
}

void popListMultiset(Multiset * multi) {
	if (!multi->done)
		chooseCoords(multi);
	if (!multi->done)
		popMultiplexers(multi);
}

void seekMultiset(Multiset * multi, const char * chrom, int start, int finish) {
	int i;
	for (i=0; i<multi->count; i++)
		seekMultiplexer(multi->multis[i], chrom, start, finish);
	multi->chrom = NULL;
	multi->start = 0;
	popMultiset(multi);
}

Multiset * newMultiset(Multiplexer ** multis, int count) {
	Multiset * new = (Multiset *) calloc (1, sizeof(Multiset));
	new->pop = popListMultiset;
	new->count = count;
	new->multis = multis;
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double **) calloc(count, sizeof(double));
	int i;
	for (i = 0; i < count; i++)
		new->values[i] = multis[i]->values;
	popMultiset(new);
	return new;
}
