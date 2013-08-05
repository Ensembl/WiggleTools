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

#include "multiplexer.h"

typedef struct combiner_st {
	char * chrom;
	int start;
	int finish;
	double * values;
	int count;
	bool *inplay;
	Multiplexer ** multis;
	bool done;
} Combiner;

void popCombiner(Combiner * combi) {
	int i, first = 0;
	int lastFinish = combi->finish;
	char * lastChrom = combi->chrom;
	bool found = false;
	int finishedIters = 0;
	int j;

	if (combi->done)
		return;

	combi->chrom = NULL;

	// Choose chromosome:
	for (i = 0; i < combi->count; i++) {
		if (combi->multis[i]->done) {
			combi->inplay[i] = false;
			finishedIters++;
		} else if (!combi->chrom || strcmp(combi->multis[i]->chrom, combi->chrom) < 0) {
			combi->chrom = combi->multis[i]->chrom;
			combi->inplay[i] = true;
			first = i;
		} else if (strcmp(combi->multis[i]->chrom, combi->chrom) > 0){
			combi->inplay[i] = false;
		} else {
			combi->inplay[i] = true;
		}
	}

	// Check for deaths
	if (finishedIters == combi->count) {
		combi->done = true;
		return;
	}

	// Reset last start if you change chromosome
	if (lastChrom && strcmp(combi->chrom, lastChrom) > 0)
		lastFinish = -1;

	// Choose start:
	for (i = 0; i < combi->count; i++) {
		int start = combi->multis[i]->start;
		if (!combi->inplay[i]) {
			continue;
		} else if (i < first) {
			combi->inplay[i] = false;
		} else if (start < lastFinish) {
			if (!found || lastFinish < combi->start) {
				combi->start = lastFinish;
				first = i;
				found = true;
			}
		} else if (start < combi->start || !found) {
			combi->start = start;
			first = i;
			found = true;
		} else if (start > combi->start) {
			combi->inplay[i] = false;
		}
	}

	// Choose finish:
	for (i = 0; i < combi->count; i++) {
		if (!combi->inplay[i]) 
			continue;
		else if (i < first) {
			combi->inplay[i] = false;
		} else if (combi->multis[i]->finish < combi->finish || i == first)
			combi->finish = combi->multis[i]->finish;
	}

	// Record values, pop as appropriate:
	for (i = 0; i < combi->count; i++) {
		if (combi->inplay[i]) {
			for (j = 0; j < combi->count) {
				if (combi->multis[i]->inplay[j]) {
					combi->values[j] = combi->multis[i]->value[j];
				}
			}
			if (combi->multis[i]->finish == combi->finish) {
				pop(combi->multis[i]);
			}
		}
	}
}

Combiner * newCombiner(Multiplexer ** multis, int count) {
	Combiner * new = (Combiner *) calloc (1, sizeof(Combiner));
	new->count = count;
	new->iters = multis;
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double *) calloc(count, sizeof(double));
	popCombiner(new);
	return new;
}
