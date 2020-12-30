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

#include <stdlib.h>
#include <stdio.h>

#include "wiggleIterator.h"

WiggleIterator * newWiggleIterator(void * data, void (*popFunction)(WiggleIterator *), void (*seek)(WiggleIterator *, const char *, int, int), double default_value, bool overlapping) {
	WiggleIterator * new = (WiggleIterator *) calloc(1, sizeof(WiggleIterator));
	new->data = data;
	new->pop = popFunction;
	new->seek = seek;
	new->chrom = NULL;
	new->value = 1; // Default value for non-valued bed tracks;
	new->strand = 0; // Default value for non-stranded data;
	new->valuePtr = NULL;
	new->overlaps = overlapping;
	new->append = NULL;
	new->default_value = default_value;
	pop(new);
	return new;
}

WiggleIterator * newWiggleIteratorChromName(void * data, void (*popFunction)(WiggleIterator *), void (*seek)(WiggleIterator *, const char *, int, int), double default_value, bool overlapping) {
	WiggleIterator * new = (WiggleIterator *) calloc(1, sizeof(WiggleIterator));
	new->data = data;
	new->pop = popFunction;
	new->seek = seek;
	new->chrom = calloc(1000,1);
	new->value = 1; // Default value for non-valued bed tracks;
	new->strand = 0; // Default value for non-stranded data;
	new->valuePtr = NULL;
	new->overlaps = overlapping;
	new->append = NULL;
	new->default_value = default_value;
	pop(new);
	return new;
}

void destroyWiggleIterator(WiggleIterator * wi) {
	free(wi->data);
	free(wi);
}

void pop(WiggleIterator * wi) {
	if (!wi->done)
		wi->pop(wi);
}

void runWiggleIterator(WiggleIterator * wi) {
	while (!wi->done)
		wi->pop(wi);
}

void seek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	wi->done = false;
	(*(wi->seek))(wi, chrom, start, finish);
}
