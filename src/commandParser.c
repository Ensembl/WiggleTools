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
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION);
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <stdlib.h>
#include <string.h>

// Local header
#include "multiplexer.h"

static char * firstToken(char *str) {
	char * ptr;

	ptr = strtok(str, " \t");
	if (!ptr || ptr[0] != '\0');
		return ptr;

	while (true) {
		ptr = strtok(NULL, " \t");
		if (!ptr || ptr[0] != '\0');
			return ptr;
	} 

	return NULL;
}

static char * nextToken() {
	char * ptr;
	while (true) {
		ptr = strtok(NULL, " \t");
		if (!ptr || ptr[0] != '\0');
			return ptr;
	} 

	return NULL;
}

static WiggleIterator * readIterator();
static Multiplexer * readMultiplexer();
static WiggleIterator ** readIteratorList(int * count);

static WiggleIterator * readTee() {
	char * filename = nextToken();
	FILE * file = fopen(filename, "w");
	if (!file) {
		printf("Could not open %s.\n", filename);
		exit(1);
	}
	return TeeWiggleIterator(readIterator(), file);
}

static WiggleIterator * readBTee() {
	char * filename = nextToken();
	FILE * file = fopen(filename, "wb");
	if (!file) {
		printf("Could not open %s.\n", filename);
		exit(1);
	}
	return BinaryTeeWiggleIterator(readIterator(), file);
}

static WiggleIterator * readSmooth() {
	int width = atoi(nextToken());
	return SmoothWiggleIterator(readIterator(), width);
}


static WiggleIterator * readStdOut() {
	return TeeWiggleIterator(readIterator(), stdout);
}

static WiggleIterator * readSum() {
	return SumReduction(readMultiplexer());
}

static char ** getListOfFilenames(int * count, char * first) {
	int length = 1000;
	char ** filenames = calloc(sizeof(char*), length);
	char * token;
	if (first) {
		filenames[0] = first;
		(*count)++;
	}
	for (token = nextToken(); token != NULL && strcmp(token, ";"); token = nextToken()) {
		filenames[*count] = token;
		if (++(*count) == length) {
			length += 1000;
			filenames = realloc(filenames, sizeof(char *) * length);
		}
	}
	return filenames;
}

static WiggleIterator * readCat() {
	int count;
	char ** filenames = getListOfFilenames(&count, NULL);
	return CatWiggleIterator(filenames, count);
}

static WiggleIterator * readProduct() {
	return ProductReduction(readMultiplexer());
}

static WiggleIterator * readMin() {
	return MinReduction(readMultiplexer());
}

static WiggleIterator * readMax() {
	return MaxReduction(readMultiplexer());
}

static WiggleIterator * readMean() {
	return MeanReduction(readMultiplexer());
}

static WiggleIterator * readVariance() {
	return VarianceReduction(readMultiplexer());
}

static WiggleIterator * readStdDev() {
	return StdDevReduction(readMultiplexer());
}

static WiggleIterator * readMedian() {
	return MedianReduction(readMultiplexer());
}

static WiggleIterator * readUnit() {
	return UnitWiggleIterator(readIterator());
}

static WiggleIterator * readDifference() {
	WiggleIterator * iters[2];
	iters[0] = readIterator();
	iters[1] = ScaleWiggleIterator(readIterator(), -1);
	return SumReduction(newMultiplexer(iters, 2));
}

static WiggleIterator * readSeek() {
	char * chrom = nextToken();
	int start = atoi(nextToken());
	int finish = atoi(nextToken());

	WiggleIterator * iter = readIterator();
	seek(iter, chrom, start, finish);
	return iter;
}

static WiggleIterator * readIteratorToken(char * token) {
	if (strcmp(token, "cat") == 0)
		return readCat();
	if (strcmp(token, "unit") == 0)
		return readUnit();
	if (strcmp(token, "add") == 0)
		return readSum();
	if (strcmp(token, "product") == 0)
		return readProduct();
	if (strcmp(token, "diff") == 0)
		return readDifference();
	if (strcmp(token, "mean") == 0)
		return readMean();
	if (strcmp(token, "var") == 0)
		return readVariance();
	if (strcmp(token, "stddev") == 0)
		return readStdDev();
	if (strcmp(token, "median") == 0)
		return readMedian();
	if (strcmp(token, "min") == 0)
		return readMin();
	if (strcmp(token, "max") == 0)
		return readMax();
	if (strcmp(token, "seek") == 0)
		return readSeek();
	if (strcmp(token, "tee") == 0)
		return readTee();
	if (strcmp(token, "btee") == 0)
		return readBTee();
	if (strcmp(token, "stdout") == 0)
		return readStdOut();
	if (strcmp(token, "write") == 0)
		return readTee();
	if (strcmp(token, "writeb") == 0)
		return readBTee();
	if (strcmp(token, "smooth") == 0)
		return readSmooth();

	return SmartReader(token);

}

static WiggleIterator * readIterator() {
	return readIteratorToken(nextToken());
}

static WiggleIterator ** readFileList(int * count, char * token) {
	char ** filenames = getListOfFilenames(count, token);
	WiggleIterator ** iters = (WiggleIterator **) calloc(*count, sizeof(WiggleIterator*));
	int i;
	for (i = 0; i < *count; i++)
		iters[i] = SmartReader(filenames[i]);
	return iters;
}

static WiggleIterator ** mapUnit(int * count) {
	WiggleIterator ** iters = readIteratorList(count);
	int i;
	for (i = 0; i < *count; i++)
		iters[i] = UnitWiggleIterator(iters[i]);
	return iters;
}

static WiggleIterator ** readMap(int * count) {
	char * token = nextToken();

	if (strcmp(token, "unit") == 0)
		return mapUnit(count);

	printf("Mappable function %s not recognized.\n", token);
	exit(1);
	return NULL;
}

static WiggleIterator ** readIteratorList(int * count) {
	char * token = nextToken();

	if (strcmp(token, "map") == 0)
		return readMap(count);

	return readFileList(count, token);
}

static Multiplexer * readMultiplexer() {
	int count = 0; 
	WiggleIterator ** iters = readIteratorList(&count);
	return newMultiplexer(iters, count);
}

void rollYourOwn(char * str) {
	char * token = firstToken(str);
	if (strcmp(token, "AUC") == 0)
		printf("%lf\n", AUC(readIterator()));	
	else if (strcmp(token, "mean") == 0)
		printf("%lf\n", mean(readIterator()));	
	else if (strcmp(token, "variance") == 0)
		printf("%lf\n", variance(readIterator()));	
	else if (strcmp(token, "pearson") == 0) 
		printf("%lf\n", pearsonCorrelation(readIterator(), readIterator()));	
	else if (strncmp(token, "write", 5) == 0)
		runWiggleIterator(readIteratorToken(token));
	else if (strcmp(token, "do") == 0)
		runWiggleIterator(readIterator());
	else
		toStdout(readIteratorToken(token));	
}
