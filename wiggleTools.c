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
	puts("\tThe program takes in Wig, BigWig, BedGraph, Bed and BigBed files, which are distinguished thanks to their suffix (.wig, .bw, bg, .bed, and .bb respectively).");
	puts("\tNote that Bed and BigBed files are treated as a binary {0,1} signal that indicates the union of regions defined in the bed.");
	puts("\tAlso note that wig, bed and bg files have to be sorted (bw and bb are already sorted)");
	puts("");
	puts("Outputs:");
	puts("\tThe program outputs a bedGraph flat file in stdout.");
	puts("");
	puts("Parameters:");
	puts("\t// Unary operators");
	puts("\twiggletools unit file");
	puts("\twiggletools abs file");
	puts("\twiggletools exp file");
	puts("\twiggletools log file");
	puts("\t");
	puts("\t// Operators between a signal and a scalar");
	puts("\twiggletools scale file factor");
	puts("\twiggletools pow file exponent");
	puts("\twiggletools exp file radix");
	puts("\twiggletools log file base");
	puts("\t");
	puts("\t// Reduction operators");
	puts("\twiggletools add file1 file2 ... ");
	puts("\twiggletools mult file1 file2 ...");
	puts("\twiggletools min file1 file2 ...");
	puts("\twiggletools max file1 file2 ...");
	puts("\twiggletools mean file1 file2 ...");
	puts("\twiggletools var file1 file2 ...");
	puts("\twiggletools stddev file1 file2 ...");
	puts("\twiggletools median file1 file2 ...");
	puts("\t");
	puts("\t// Calculations");
	puts("\twiggletools AUC file");
	puts("\twiggletools pearson file1 file2");
	puts("");
	puts("\t// Other");
	puts("\twiggletools --help");
}

int main(int argc, char ** argv) {
	if (argc < 2 || strcmp(argv[1], "help") == 0) {
		printHelp();
		return 0;
	} else if (strcmp(argv[1], "add") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(SumWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "scale") == 0) {
		toStdout(ScaleWiggleIterator(WigOrBigWigReader(argv[2]), atoi(argv[3])));
	} else if (strcmp(argv[1], "mult") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(ProductWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "min") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(MinWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "max") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(MaxWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "mean") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(MeanWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "median") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(MedianWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "stddev") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(StdDevWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "var") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		toStdout(VarianceWiggleReducer(iters, count));
	} else if (strcmp(argv[1], "pow") == 0) {
		toStdout(PowerWiggleIterator(WigOrBigWigReader(argv[2]), atoi(argv[3])));
	} else if (strcmp(argv[1], "exp") == 0) {
		if (argc == 4)
			toStdout(ExpWiggleIterator(WigOrBigWigReader(argv[2]), atoi(argv[3])));
		else
			toStdout(NaturalExpWiggleIterator(WigOrBigWigReader(argv[2])));
	} else if (strcmp(argv[1], "log") == 0) {
		if (argc == 4)
			toStdout(LogWiggleIterator(WigOrBigWigReader(argv[2]), atoi(argv[3])));
		else
			toStdout(NaturalLogWiggleIterator(WigOrBigWigReader(argv[2])));
	} else if (strcmp(argv[1], "unit") == 0) 
		toStdout(UnitWiggleIterator(WigOrBigWigReader(argv[2])));
	else if (strcmp(argv[1], "abs") == 0) 
		toStdout(AbsWiggleIterator(WigOrBigWigReader(argv[2])));
	else if (strcmp(argv[1], "--help") == 0) 
		printHelp();
	else if (strcmp(argv[1], "pearson") == 0)
		printf("%f\n", pearsonCorrelation(WigOrBigWigReader(argv[2]), WigOrBigWigReader(argv[3])));
	else if (strcmp(argv[1], "AUC") == 0) 
		printf("%f\n", AUC(WigOrBigWigReader(argv[2])));
	else if (strcmp(argv[1], "stream") == 0) {
		int i;
		int count = argc - 2;
		WiggleIterator ** iters = (WiggleIterator **) calloc(count, sizeof(WiggleIterator*));
		for (i = 0; i < count; i++)
			iters[i] = WigOrBigWigReader(argv[i + 2]);
		streamMultiplexer(stdout, newMultiplexer(iters, count));
	} else if (strcmp(argv[1], "catch") == 0) {
		streamMultiplexer(stdout, newStreamingMultiplexer(stdin));
	} else {
		printf("Unrecognized keyword: %s\n", argv[1]);
		puts("");
		printHelp();
		return 1;
	}

	return 0;
}

