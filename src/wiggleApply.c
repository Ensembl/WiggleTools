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

// Local header
#include "wiggleTools.h"
#include "wiggleIterators.h"

typedef struct applyWiggleIteratorTarget_st {
	WiggleIterator * regions;
	double (*statistic)(WiggleIterator *);
	WiggleIterator * data;
	struct applyWiggleIteratorTarget_st * next;
} ApplyWiggleIteratorTarget;

typedef struct applyWiggleIteratorData_st {
	WiggleIterator * regions;
	double (*statistic)(WiggleIterator *);
	WiggleIterator * data;
} ApplyWiggleIteratorData;

void ApplyWiggleIteratorPop(WiggleIterator * wi) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) wi->data;
	wi->chrom = data->regions->chrom;
	wi->start = data->regions->start;
	wi->finish = data->regions->finish;
	seek(wi->data, wi->chrom, wi->start, wi->finish);
	wi->value = (*(data->statistic))(data->data);
}

void ApplyWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	seek(wi->regions, chrom, start, finish);
}

WiggleIterator * apply(WiggleIterator * regions, double (*statistic)(WiggleIterator *), WiggleIterator * dataset) {
	ApplyWiggleIteratorData * data = (ApplyWiggleIteratorData *) calloc(1, sizeof(ApplyWiggleIteratorData));
	data->regions = regions;
	data->statistic = statistic;
	data->data = dataset;
	return newWiggleIterator(data, &ApplyWiggleIteratorPop, &ApplyWiggleIteratorSeek);
}
