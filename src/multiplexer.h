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

#ifndef WIGGLE_MULTIPLEXER_H_
#define WIGGLE_MULTIPLEXER_H_

#include "wiggleIterator.h"
#include "fib.h"

struct multiplexer_st {
	char * chrom;
	int start;
	int finish;
	double * values;
	double * default_values;
	int count, inplay_count;
	bool *inplay;
	WiggleIterator ** iters;
	bool done;
	bool strict;
	void (*pop)(Multiplexer *);
	void (*seek)(Multiplexer *, const char *, int, int);
	FibHeap * starts, *finishes;
	void * data;
};

void popMultiplexer(Multiplexer * multi);
void seekMultiplexer(Multiplexer * multi, const char * chrom, int start, int finish);
void runMultiplexer(Multiplexer * multi);
Multiplexer * newCoreMultiplexer(void * data, int count, void (*pop)(Multiplexer *), void (*seek)(Multiplexer *, const char *, int, int));

#endif
