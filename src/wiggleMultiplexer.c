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
			*inplayPtr = true;
			if (multi->start > clipping && start < multi->start) {
				if (start < clipping)
					multi->start = clipping;
				else
					multi->start = start;
				multi->finish = (*wiPtr)->finish;
				first = i;
			} else if (start > multi->start) {
				*inplayPtr = false;
			} else if ((*wiPtr)->finish < multi->finish) {
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

static void readValues(Multiplexer * multi) {
	int i; 
	bool * inplayPtr = multi->inplay;
	WiggleIterator ** wiPtr = multi->iters;
	double * valuePtr = multi->values;

	for (i=0; i < multi->count; i++) {
		if (*inplayPtr) {
			*valuePtr = (*wiPtr)->value;
			if ((*wiPtr)->finish == multi->finish) {
				pop(*wiPtr);
			}
		}
		valuePtr++;
		inplayPtr++;
		wiPtr++;
	}
}

void popListMultiplexer(Multiplexer * multi) {
	if (multi->done)
		return;

	chooseCoords(multi);

	if (multi->done) 
		return;

	readValues(multi);

}

void seekMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish) {
	int i;
	for (i=0; i<multi->count; i++)
		seek(multi->iters[i], chrom, start, finish);
	multi->chrom = NULL;
	multi->start = 0;
	popMultiplexer(multi);
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
