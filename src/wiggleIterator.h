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

#ifndef _WIGGLETOOLS_PRIV_
#define _WIGGLETOOLS_PRIV_

#include <stdio.h>
#include "wiggletools.h"

struct wiggleIterator_st {
	char * chrom;
	int start;
	int finish;
	double value;
	void * valuePtr;
	bool done;
	int strand;
	void * data;
	void (*pop)(WiggleIterator *);
	void (*seek)(WiggleIterator *, const char *, int, int);
	bool overlaps;
	double default_value;
};

WiggleIterator * newWiggleIterator(void * data, void (*pop)(WiggleIterator *), void (*seek)(WiggleIterator *, const char *, int, int));
void pop(WiggleIterator *);
WiggleIterator * CompressionWiggleIterator(WiggleIterator *);

#endif
