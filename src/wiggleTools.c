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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
// Local header
#include "wiggleTools.h"

static void printHelp() {
	puts("Inputs:");
	puts("\tThe program takes in Wig, BigWig, BedGraph, Bed, BigBed and Bam files, which are distinguished thanks to their suffix (.wig, .bw, bg, .bed, .bb, .bam respectively).");
	puts("\tNote that Bed and BigBed files are treated as a binary {0,1} signal that indicates the union of regions defined in the bed.");
	puts("\tAlso note that .wig, .bed and .bg files have to be sorted by coordinate (.bw, .bb and .bam files are already internally sorted)");
	puts("\tFinally, the program assumes that every .bam file has an ancillary .bai index file in the same directory");
	puts("");
	puts("Outputs:");
	puts("\tThe program outputs a bedGraph flat file in stdout.");
	puts("");
	puts("Parameters:");
	puts("\t// Unary operators");
	puts("\twiggletools unit [options] file");
	puts("\twiggletools abs [options] file");
	puts("\twiggletools exp [options] file");
	puts("\twiggletools log [options] file");
	puts("\t");
	puts("\t// Operators between a signal and a scalar");
	puts("\twiggletools scale [options] file factor");
	puts("\twiggletools pow [options] file exponent");
	puts("\twiggletools exp [options] file radix");
	puts("\twiggletools log [options] file base");
	puts("\t");
	puts("\t// Reduction operators");
	puts("\twiggletools add [options] file1 file2 ... ");
	puts("\twiggletools mult [options] file1 file2 ...");
	puts("\twiggletools min [options] file1 file2 ...");
	puts("\twiggletools max [options] file1 file2 ...");
	puts("\twiggletools mean [options] file1 file2 ...");
	puts("\twiggletools var [options] file1 file2 ...");
	puts("\twiggletools stddev [options] file1 file2 ...");
	puts("\twiggletools median [options] file1 file2 ...");
	puts("\t");
	puts("\t// Calculations");
	puts("\twiggletools AUC [options] file");
	puts("\twiggletools pearson [options] file1 file2");
	puts("");
	puts("\t// Applied statistics");
	puts("\twiggletools coverage regions data");
	puts("\t// Other");
	puts("\twiggletools --help");
}

//static int stripOptions(int argc, char ** argv) {
//	int i;
//
//	for (i = 0; (i < argc) && (argv[i][0] == '-'); i++)
//		parseOption(argv[i]);
//
//	return i;
//}

static WiggleIterator ** SmartReaders(char ** filenames, int count) {
	int i;
	WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
	for (i = 0; i < count; i++)
		iters[i] = SmartReader(filenames[i]);
	return iters;
}

int main(int argc, char ** argv) {
	if (argc < 2 || strcmp(argv[1], "help") == 0) {
		printHelp();
		return 0;
	} else if (strcmp(argv[1], "add") == 0) {
		toStdout(SumWiggleReducer(SmartReaders(argv + 2, argc - 2), argc - 2));
	} else if (strcmp(argv[1], "scale") == 0) {
		toStdout(ScaleWiggleIterator(SmartReader(argv[2]), atoi(argv[3])));
	} else if (strcmp(argv[1], "mult") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(ProductWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "min") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(MinWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "max") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(MaxWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "mean") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(MeanWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "median") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(MedianWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "stddev") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(StdDevWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "var") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		toStdout(VarianceWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "pow") == 0) {
		toStdout(PowerWiggleIterator(SmartReader(argv[2]), atoi(argv[3])));
	} else if (strcmp(argv[1], "exp") == 0) {
		if (argc == 4)
			toStdout(ExpWiggleIterator(SmartReader(argv[2]), atoi(argv[3])));
		else
			toStdout(NaturalExpWiggleIterator(SmartReader(argv[2])));
	} else if (strcmp(argv[1], "log") == 0) {
		if (argc == 4)
			toStdout(LogWiggleIterator(SmartReader(argv[2]), atoi(argv[3])));
		else
			toStdout(NaturalLogWiggleIterator(SmartReader(argv[2])));
	} else if (strcmp(argv[1], "unit") == 0) 
		toStdout(UnitWiggleIterator(SmartReader(argv[2])));
	else if (strcmp(argv[1], "abs") == 0) 
		toStdout(AbsWiggleIterator(SmartReader(argv[2])));
	else if (strcmp(argv[1], "--help") == 0) 
		printHelp();
	else if (strcmp(argv[1], "pearson") == 0)
		printf("%f\n", pearsonCorrelation(SmartReader(argv[2]), SmartReader(argv[3])));
	else if (strcmp(argv[1], "AUC") == 0) 
		printf("%f\n", AUC(SmartReader(argv[2])));
	else if (strcmp(argv[1], "coverage") == 0) 
		toStdout(apply(SmartReader(argv[2]), AUC, SmartReader(argv[3])));
	else if (strcmp(argv[1], "stream") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = SmartReader(argv[i + 2]);
		streamMultiplexer(stdout, newMultiplexer(iters, count));
	} else if (strcmp(argv[1], "catch") == 0) {
		streamMultiplexer(stdout, newStreamingMultiplexer(stdin));
	} else if (strcmp(argv[1], "seek") == 0) {
		WiggleIterator * wi = SmartReader(argv[2]);
		seek(wi, argv[3], atoi(argv[4]), atoi(argv[5]));
		toStdout(wi);
	} else {
		printf("Unrecognized keyword: %s\n", argv[1]);
		puts("");
		printHelp();
		return 1;
	}

	return 0;
}

