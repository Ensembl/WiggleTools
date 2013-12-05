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

#include <stdlib.h>
#include <stdio.h>

#include "wiggleIterator.h"

WiggleIterator * newWiggleIterator(void * data, void (*popFunction)(WiggleIterator *), void (*seek)(WiggleIterator *, const char *, int, int)) {
	WiggleIterator * new = (WiggleIterator *) calloc(1, sizeof(WiggleIterator));
	new->data = data;
	new->pop = popFunction;
	new->seek = seek;
	new->chrom = calloc(1000,1);
	new->value = 1; // Default value for non-valued bed tracks;
	new->strand = 0; // Default value for non-stranded data;
	new->valuePtr = NULL;
	new->overlaps = false;
	pop(new);
	return new;
}

void destroyWiggleIterator(WiggleIterator * wi) {
	free(wi->data);
	free(wi->chrom);
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
	(*(wi->seek))(wi, chrom, start, finish);
}
