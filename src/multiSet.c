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
