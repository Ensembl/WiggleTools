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
#include <string.h>

// Local header
#include "multiplexer.h"

void printHelp() {

puts("WiggleTools");
puts("");
puts("Copyright EMBL-EBI, 2013.");
puts("Development contact: Daniel Zerbino zerbino@ebi.ac.uk");
puts("");
puts("This library parses wiggle files and executes various operations on them streaming through lazy evaluators.");
puts("");
puts("Inputs:");
puts("\tThe program takes in Wig, BigWig, BedGraph, Bed, BigBed and Bam files, which are distinguished thanks to their suffix (.wig, .bw, .bg, .bed, .bb, .bam respectively).");
puts("\tNote that wiggletools assumes that every bam file has an index .bai file next to it.");
puts("");
puts("Outputs:");
puts("\tThe program outputs a wiggle file in stdout unless the output is squashed");
puts("");
puts("Command line:");
puts("\twiggletools --help");
puts("\twiggletools ' program '");
puts("");
puts("Program grammar:");
puts("\tprogram = (iterator) | do (iterator) | (statistic) | (extraction)");
puts("\tstatistic = AUC (iterator) | mean (iterator) | variance (iterator) | pearson (iterator) (iterator) | isZero (iterator)");
puts("\textraction = profile (iterator) (iterator) | profiles (iterator) (iterator)");
puts("\t\t| apply (out_filename) (statistic) (bed_file) (iterator)");
puts("\titerator = (filename) | (unary_operator) (iterator) | (binary_operator) (iterator) (iterator) | (reducer) (multiplex) | (setComparison) (multiplex) (multiplex)");
puts("\tunary_operator = unit | stdout | write (filename.wig) | smooth (int) | exp | ln | log (double) | pow (double)");
puts("\tbinary_operator = diff | ratio");
puts("\tmultiplex = (filename_list) | map (unary_operator) (multiplex)");
puts("\treducer = cat | sum | product | mean | var | stddev | median | min | max");
puts("\tsetComparison = ttest | wilcoxon");
puts("\tfilename_list = (filename) : | (filename) (filename_list)");
puts("\tfilename = *.wig | *.bw | *.bed | *.bb | *.bg | *.bam");

}

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

static WiggleIterator * readApply() {
	char * outfilename = nextToken();
	char * operation = nextToken();
	char * infilename = nextToken();

	FILE * infile = fopen(infilename, "r");
	if (!infile) {
		printf("Could not open %s.\n", infilename);
		exit(1);
	}

	FILE * outfile;
	if (strcmp(outfilename, "-"))
		outfile = fopen(outfilename, "w");
	else
		outfile = stdout;

	if (!outfile) {
		printf("Could not open %s.\n", outfilename);
		exit(1);
	}

	double (*function)(WiggleIterator *);

	if (strcmp(operation, "AUC") == 0)
		function = &AUC;
	else if (strcmp(operation, "mean") == 0)
		function = &mean;
	else if (strcmp(operation, "variance") == 0)
		function = &variance;
	else {
		printf("Name of function to be applied unrecognized: %s\n", operation);
		exit(1);
	}

	return PasteWiggleIterator(ApplyWiggleIterator(SmartReader(infilename), function, readIterator()), infile, outfile);
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

static WiggleIterator * readPow() {
	double base = atof(nextToken());
	return PowerWiggleIterator(readIterator(), base);
}

static WiggleIterator * readGt() {
	double cutoff = atof(nextToken());
	return HighPassFilterWiggleIterator(readIterator(), cutoff);
}

static WiggleIterator * readScale() {
	double scalar = atof(nextToken());
	return ScaleWiggleIterator(readIterator(), scalar);
}

static WiggleIterator * readExp() {
	return NaturalExpWiggleIterator(readIterator());
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

static WiggleIterator * readNaturalLog() {
	return NaturalLogWiggleIterator(readIterator());
}

static WiggleIterator * readLog() {
	double base = atof(needNextToken());
	return LogWiggleIterator(readIterator(), base);
}

static WiggleIterator * readTTest() {
	Multiplexer ** multis = calloc(2, sizeof(Multiplexer *));
	multis[0] = readMultiplexer();
	multis[1] = readMultiplexer();
	return TTestReduction(newMultiset(multis, 2));
}

static WiggleIterator * readMWUTest() {
	Multiplexer ** multis = calloc(2, sizeof(Multiplexer *));
	multis[0] = readMultiplexer();
	multis[1] = readMultiplexer();
	return MWUReduction(newMultiset(multis, 2));
}

static WiggleIterator * readIteratorToken(char * token) {
	if (strcmp(token, "cat") == 0)
		return readCat();
	if (strcmp(token, "scale") == 0)
		return readScale();
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
	if (strcmp(token, "exp") == 0)
		return readExp();
	if (strcmp(token, "ln") == 0)
		return readNaturalLog();
	if (strcmp(token, "log") == 0)
		return readLog();
	if (strcmp(token, "pow") == 0)
		return readPow();
	if (strcmp(token, "gt") == 0)
		return readGt();
	if (strcmp(token, "ttest") == 0)
		return readTTest();
	if (strcmp(token, "wilcoxon") == 0)
		return readMWUTest();

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

static void readProfile() {
	int width = atoi(nextToken());
	char * regions = nextToken();
	WiggleIterator * wig = readIterator();
	double * profile = profileSum(SmartReader(regions), wig, width);
	int i;
	for (i = 0; i < width; i++)
		printf("%f\n", profile[i]);
	free(profile);
}

static void readProfiles() {
	int width = atoi(nextToken());
	double * profile = calloc(width, sizeof(double));
	WiggleIterator * regions = readIterator();
	WiggleIterator * wig = readIterator();

	for (; !regions->done; pop(regions)) {
		int i;
		regionProfile(regions, wig, width, profile, width/2);
		printf("%f", profile[0]);
		profile[0] = 0;

		for (i = 1; i < width; i++) {
			printf("\t%f", profile[i]);
			profile[i] = 0;
		}

		printf("\n");
	}

	free(profile);
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
	else if (strcmp(token, "isZero") == 0)
		isZero(readIterator());	
	else if (strcmp(token, "apply") == 0)
		runWiggleIterator(readApply());
	else if (strcmp(token, "profile") == 0)
		readProfile();
	else if (strcmp(token, "profiles") == 0)
		readProfiles();
	else
		toStdout(readIteratorToken(token));	
}
