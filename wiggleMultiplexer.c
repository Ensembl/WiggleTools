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
//

#include <stdlib.h>
#include <string.h>

#include "wiggleMultiplexer.h"

void popMultiplexer(Multiplexer * multi) {
	(*(multi->pop))(multi);
}

void popListMultiplexer(Multiplexer * multi) {
	int i, first = 0;
	int lastFinish = multi->finish;
	char * lastChrom = multi->chrom;
	bool found = false;
	int finishedIters = 0;

	if (multi->done)
		return;

	multi->chrom = NULL;

	// Choose chromosome:
	for (i = 0; i < multi->count; i++) {
		if (multi->iters[i]->done) {
			multi->inplay[i] = false;
			finishedIters++;
		} else if (!multi->chrom || strcmp(multi->iters[i]->chrom, multi->chrom) < 0) {
			multi->chrom = multi->iters[i]->chrom;
			multi->inplay[i] = true;
			first = i;
		} else if (strcmp(multi->iters[i]->chrom, multi->chrom) > 0){
			multi->inplay[i] = false;
		} else {
			multi->inplay[i] = true;
		}
	}

	// Check for deaths
	if (finishedIters == multi->count) {
		multi->done = true;
		return;
	}

	// Reset last start if you change chromosome
	if (lastChrom && strcmp(multi->chrom, lastChrom) > 0)
		lastFinish = -1;

	// Choose start:
	for (i = 0; i < multi->count; i++) {
		int start = multi->iters[i]->start;
		if (!multi->inplay[i]) {
			continue;
		} else if (i < first) {
			multi->inplay[i] = false;
		} else if (start < lastFinish) {
			if (!found || lastFinish < multi->start) {
				multi->start = lastFinish;
				first = i;
				found = true;
			}
		} else if (start < multi->start || !found) {
			multi->start = start;
			first = i;
			found = true;
		} else if (start > multi->start) {
			multi->inplay[i] = false;
		}
	}

	// Choose finish:
	for (i = 0; i < multi->count; i++) {
		if (!multi->inplay[i]) 
			continue;
		else if (i < first) {
			multi->inplay[i] = false;
		} else if (multi->iters[i]->finish < multi->finish || i == first)
			multi->finish = multi->iters[i]->finish;
	}

	// Record values, pop as appropriate:
	for (i = 0; i < multi->count; i++) {
		if (multi->inplay[i]) {
			multi->values[i] = multi->iters[i]->value;
			if (multi->iters[i]->finish == multi->finish) {
				pop(multi->iters[i]);
			}
		}
	}

}

void seekMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	int i;
	for (i=0; i<multi->count; i++)
		seek(multi->iters[i], chrom, start, finish);
}

Multiplexer * newMultiplexer(WiggleIterator ** iters, int count) {
	Multiplexer * new = (Multiplexer *) calloc (1, sizeof(Multiplexer));
	new->pop = popListMultiplexer;
	new->count = count;
	new->iters = iters;
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double *) calloc(count, sizeof(double));
	popMultiplexer(new);
	return new;
}

Multiplexer * newIteratorMultiplexer(WiggleIterator * iter, int index, int count) {
	int i;

	Multiplexer * new = (Multiplexer *) calloc (1, sizeof(Multiplexer));
	new->pop = popListMultiplexer;
	new->count = count;
	new->iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator *));
	for (i = 0; i< count; i++) {
		if (i != count)
			new->iters[i] = NullWiggleIterator();
		else
			new->iters[i] = iter;
	}
	new->inplay = (bool *) calloc(count, sizeof(bool));
	new->values = (double *) calloc(count, sizeof(double));
	popMultiplexer(new);
	return new;
}
