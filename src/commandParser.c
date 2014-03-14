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
puts("\twiggletools program");
puts("");
puts("Program grammar:");
puts("\tprogram = (iterator) | do (iterator) | (statistic) | (extraction)");
puts("\tstatistic = AUC (output) (iterator) | meanI (output) (iterator) | varI (output) (iterator) | pearson (output) (iterator) (iterator) | isZero (iterator)");
puts("\toutput = filename | -");
puts("\textraction = profile (output) (int) (iterator) (iterator) | profiles (output) (int) (iterator) (iterator) | histogram (output) (width)");
puts("\t\t| apply_paste (out_filename) (statistic) (bed_file) (iterator)");
puts("\titerator = (filename) | (unary_operator) (iterator) | (binary_operator) (iterator) (iterator) | (reducer) (multiplex) | (setComparison) (multiplex) (multiplex)");
puts("\tunary_operator = unit | write (output) | write_bg (ouput) | smooth (int) | exp | ln | log (float) | pow (float) | offset (float) | scale (float) | gt (float)");
puts("\tbinary_operator = diff | ratio | overlaps | apply (statistic)");
puts("\treducer = cat | sum | product | mean | var | stddev | entropy | CV | median | min | max");
puts("\titerator_list = (iterator) : | (iterator) (iterator_list)");
puts("\tmultiplex = (iterator_list) | map (unary_operator) (multiplex)");
puts("\tsetComparison = ttest | wilcoxon");
puts("\tfilename = *.wig | *.bw | *.bed | *.bb | *.bg | *.bam");

}

static char * nextToken(int argc, char ** argv) {
	static int count;
	static char ** ptr;
	static int index = 0;
	if (argv) {
		ptr = argv;
		count = argc;
	}
	if (index == count)
		return NULL;
	else
		return ptr[index++];
}

static char * needNextToken() {
	char * token = nextToken(0,0);
	if (token) {
		return token;
	} else {
		fprintf(stderr, "wiggletools: Unexpected end of command line\n");
		exit(1);
	}
}

static WiggleIterator * readIteratorToken(char * token);

static WiggleIterator * readIterator() {
	return readIteratorToken(needNextToken());
}

static WiggleIterator ** readFileList(int * count, char * firstToken) {
	size_t buffer_size = 8;
	char * token;
	int i =0;
	WiggleIterator ** iters = (WiggleIterator **) calloc(buffer_size, sizeof(WiggleIterator*));

	for (token = firstToken; token != NULL && strcmp(token, ":"); token = nextToken(0,0)) {
		if (i == buffer_size) {
			buffer_size *= 2;
			iters = (WiggleIterator **) realloc(iters, buffer_size * sizeof(WiggleIterator*));
		}
		iters[i++] = readIteratorToken(token);
	}
	*count = i;
	return iters;
}

static WiggleIterator ** readIteratorList(int * count);

static WiggleIterator ** readMappedIteratorList(int * count) {
	char * token = needNextToken();
	WiggleIterator ** iters;
	int i;

	if (strcmp(token, "unit") == 0) {
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = UnitWiggleIterator(iters[i]);
	} else if (strcmp(token, "smooth") == 0) {
		int width = atoi(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = SmoothWiggleIterator(iters[i], width);
	} else if (strcmp(token, "exp") == 0) {
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = NaturalExpWiggleIterator(iters[i]);
	} else if (strcmp(token, "ln") == 0) {
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = NaturalLogWiggleIterator(iters[i]);
	} else if (strcmp(token, "log") == 0) {
		double base = atof(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = LogWiggleIterator(iters[i], base);
	} else if (strcmp(token, "pow") == 0) {
		double base = atof(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = PowerWiggleIterator(iters[i], base);
	} else if (strcmp(token, "scale") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = ScaleWiggleIterator(iters[i], scalar);
	} else if (strcmp(token, "offset") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = ShiftWiggleIterator(iters[i], scalar);
	} else if (strcmp(token, "gt") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count);
		for (i = 0; i < *count; i++)
			iters[i] = HighPassFilterWiggleIterator(iters[i], scalar);
	} else {
		fprintf(stderr, "Unary function unkown: %s\n", token);
		exit(1);
	}

	return iters;
}

static WiggleIterator ** readIteratorList(int * count) {
	char * token = needNextToken();

	if (strcmp(token, "map"))
		return readFileList(count, token);
	else 
		return readMappedIteratorList(count);
}

static Multiplexer * readMultiplexer() {
	int count = 0; 
	WiggleIterator ** iters = readIteratorList(&count);
	return newMultiplexer(iters, count);
}

static FILE * readOutputFilename() {
	char * filename = needNextToken();
	if (strcmp(filename, "-")) {
		FILE * file = fopen(filename, "w");
		if (!file) {
			fprintf(stderr, "Could not open output file %s.\n", filename);
			exit(1);
		}
		return file;
	} else 
		return stdout;
}

static WiggleIterator * readTee() {
	FILE * file = readOutputFilename();
	return TeeWiggleIterator(readIterator(), file, false);
}

static WiggleIterator * readBGTee() {
	FILE * file = readOutputFilename();
	return TeeWiggleIterator(readIterator(), file, true);
}

static WiggleIterator * readLastIteratorToken(char * token) {
	WiggleIterator * iter = readIteratorToken(token);
	char * remainder = nextToken(0,0);
	if (remainder) {
		fprintf(stderr, "Trailing tokens: the last tokens in your command were not read, check your syntax:\n...");
		while (remainder) {
			fprintf(stderr, " %s", remainder);
			remainder = nextToken(0,0);
		}
		fprintf(stderr, "\n");
		exit(1);
	}
	return iter;
}

static WiggleIterator * readLastIterator() {
	return readLastIteratorToken(needNextToken());
}

static WiggleIterator * readSmooth() {
	int width = atoi(needNextToken());
	return SmoothWiggleIterator(readIterator(), width);
}

static WiggleIterator * readPow() {
	double base = atof(needNextToken());
	return PowerWiggleIterator(readIterator(), base);
}

static WiggleIterator * readOverlap() {
	WiggleIterator * source = readIterator();
	WiggleIterator * mask = readIterator();
	return OverlapWiggleIterator(source, mask);
}

static WiggleIterator * readGt() {
	double cutoff = atof(needNextToken());
	return HighPassFilterWiggleIterator(readIterator(), cutoff);
}

static WiggleIterator * readScale() {
	double scalar = atof(needNextToken());
	return ScaleWiggleIterator(readIterator(), scalar);
}

static WiggleIterator * readShift() {
	double scalar = atof(needNextToken());
	return ShiftWiggleIterator(readIterator(), scalar);
}

static WiggleIterator * readExp() {
	return NaturalExpWiggleIterator(readIterator());
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
	for (token = needNextToken(); token != NULL && strcmp(token, ":"); token = nextToken(0,0)) {
		filenames[*count] = token;
		if (++(*count) == length) {
			length += 1000;
			filenames = realloc(filenames, sizeof(char *) * length);
		}
	}
	return filenames;
}

static WiggleIterator * readCat() {
	int count = 0;
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

static WiggleIterator * readEntropy() {
	return EntropyReduction(readMultiplexer());
}

static WiggleIterator * readCV() {
	return CVReduction(readMultiplexer());
}

static WiggleIterator * readMedian() {
	return MedianReduction(readMultiplexer());
}

static WiggleIterator * readUnit() {
	return UnitWiggleIterator(readIterator());
}

static WiggleIterator * readDifference() {
	WiggleIterator ** iters = calloc(2, sizeof(WiggleIterator *));
	iters[0] = readIterator();
	iters[1] = ScaleWiggleIterator(readIterator(), -1);
	return SumReduction(newMultiplexer(iters, 2));
}

static WiggleIterator * readRatio() {
	WiggleIterator ** iters = calloc(2, sizeof(WiggleIterator *));
	iters[0] = readIterator();
	iters[1] = PowerWiggleIterator(readIterator(), -1);
	return ProductReduction(newMultiplexer(iters, 2));
}

static WiggleIterator * readSeek() {
	char * chrom = needNextToken();
	int start = atoi(needNextToken());
	int finish = atoi(needNextToken());

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

static WiggleIterator * readApply() {
	char * operation = needNextToken();
	double (*function)(WiggleIterator *);

	if (strcmp(operation, "AUC") == 0)
		function = &AUC;
	else if (strcmp(operation, "mean") == 0)
		function = &mean;
	else if (strcmp(operation, "variance") == 0)
		function = &variance;
	else {
		fprintf(stderr, "Name of function to be applied unrecognized: %s\n", operation);
		exit(1);
	}

	WiggleIterator * regions = readIterator();
	WiggleIterator * data = readIterator();

	return ApplyWiggleIterator(regions, function, data);
}


static WiggleIterator * readIteratorToken(char * token) {
	if (strcmp(token, "cat") == 0)
		return readCat();
	if (strcmp(token, "scale") == 0)
		return readScale();
	if (strcmp(token, "offset") == 0)
		return readShift();
	if (strcmp(token, "unit") == 0)
		return readUnit();
	if (strcmp(token, "sum") == 0)
		return readSum();
	if (strcmp(token, "mult") == 0)
		return readProduct();
	if (strcmp(token, "diff") == 0)
		return readDifference();
	if (strcmp(token, "ratio") == 0)
		return readRatio();
	if (strcmp(token, "apply") == 0)
		return readApply();
	if (strcmp(token, "mean") == 0)
		return readMean();
	if (strcmp(token, "var") == 0)
		return readVariance();
	if (strcmp(token, "stddev") == 0)
		return readStdDev();
	if (strcmp(token, "entropy") == 0)
		return readEntropy();
	if (strcmp(token, "CV") == 0)
		return readCV();
	if (strcmp(token, "median") == 0)
		return readMedian();
	if (strcmp(token, "min") == 0)
		return readMin();
	if (strcmp(token, "max") == 0)
		return readMax();
	if (strcmp(token, "seek") == 0)
		return readSeek();
	if (strcmp(token, "write") == 0)
		return readTee();
	if (strcmp(token, "write_bg") == 0)
		return readBGTee();
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
	if (strcmp(token, "overlaps") == 0)
		return readOverlap();
	if (strcmp(token, "ttest") == 0)
		return readTTest();
	if (strcmp(token, "wilcoxon") == 0)
		return readMWUTest();

	return SmartReader(token);

}

static void readProfile() {
	FILE * file = readOutputFilename();

	int width = atoi(needNextToken());
	WiggleIterator * regions = readIterator();
	WiggleIterator * wig = readLastIterator();
	WiggleIterator * profiles = ProfileWiggleIterator(regions, width, wig);
	double * profile = calloc(width, sizeof(double));

	for (; !profiles->done; pop(profiles))
		addProfile(profile, (double *) profiles->valuePtr, width);

	int i;
	for (i = 0; i < width; i++)
		fprintf(file, "%i\t%lf\n", i, profile[i]);

	free(profile);
	fclose(file);
}

static void fprintfProfile(FILE * file, double * profile, int width) {
	int i;

	fprintf(file, "%f", profile[0]);
	for (i = 1; i < width; i++)
		fprintf(file, "\t%f", profile[i]);
	fprintf(file, "\n");
}

static void readProfiles() {
	FILE * file = readOutputFilename();

	int width = atoi(needNextToken());
	WiggleIterator * regions = readIterator();
	WiggleIterator * wig = readLastIterator();
	WiggleIterator * profiles;

	for (profiles = ProfileWiggleIterator(regions, width, wig); !profiles->done; pop(profiles)) {
		fprintf(file, "%s\t%i\t%i\t", profiles->chrom, profiles->start, profiles->finish);
		fprintfProfile(file, (double *) profiles->valuePtr, width);
	}

	fclose(file);
}

static void readAUC() {
	FILE * file = readOutputFilename();
	fprintf(file, "%lf\n", AUC(readLastIterator()));	
	fclose(file);
}
 
static void readHistogram() {
	FILE * file = readOutputFilename();
	int width = atoi(needNextToken());
	Histogram * hist = histogram(readLastIterator(), width);	
	print_histogram(hist, file);
	fclose(file);
}

static void readMeanIntegrated() {
	FILE * file = readOutputFilename();
	fprintf(file, "%lf\n", mean(readLastIterator()));	
	fclose(file);
}

static void readVarianceIntegrated() {
	FILE * file = readOutputFilename();
	fprintf(file, "%lf\n", variance(readLastIterator()));	
	fclose(file);
}

static void readPearson() {
	FILE * file = readOutputFilename();
	WiggleIterator * iter1 = readIterator();
	WiggleIterator * iter2 = readLastIterator();
	fprintf(file, "%lf\n", pearsonCorrelation(iter1, iter2));
	fclose(file);
}

static WiggleIterator * readApplyPaste() {
       FILE * outfile = readOutputFilename();
       char * operation = needNextToken();
       char * infilename = needNextToken();

       FILE * infile = fopen(infilename, "r");
       if (!infile) {
               fprintf(stderr, "Could not open %s.\n", infilename);
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
               fprintf(stderr, "Name of function to be applied unrecognized: %s\n", operation);
               exit(1);
       }

       return PasteWiggleIterator(ApplyWiggleIterator(SmartReader(infilename), function, readLastIterator()), infile, outfile);
}

void rollYourOwn(int argc, char ** argv) {
	char * token = nextToken(argc, argv);
	if (strcmp(token, "AUC") == 0)
		readAUC();
	else if (strcmp(token, "histogram") == 0)
		readHistogram();
	else if (strcmp(token, "meanI") == 0)
		readMeanIntegrated();
	else if (strcmp(token, "varI") == 0)
		readVarianceIntegrated();
	else if (strcmp(token, "pearson") == 0) 
		readPearson();
	else if (strncmp(token, "write", 5) == 0)
		runWiggleIterator(readLastIteratorToken(token));
	else if (strcmp(token, "do") == 0)
		runWiggleIterator(readLastIterator());
	else if (strcmp(token, "isZero") == 0)
		isZero(readLastIterator());	
	else if (strcmp(token, "apply_paste") == 0)
		runWiggleIterator(readApplyPaste());
	else if (strcmp(token, "profile") == 0)
		readProfile();
	else if (strcmp(token, "profiles") == 0)
		readProfiles();
	else
		toStdout(readLastIteratorToken(token), false);	
}
