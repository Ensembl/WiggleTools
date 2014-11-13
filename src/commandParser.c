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

bool holdFire = false;

void printHelp() {

puts("WiggleTools");
puts("");
puts("Copyright EMBL-EBI, 2013.");
puts("Development contact: Daniel Zerbino zerbino@ebi.ac.uk");
puts("");
puts("This library parses wiggle files and executes various operations on them streaming through lazy evaluators.");
puts("");
puts("Inputs:");
puts("\tThe program takes in Wig, BigWig, BedGraph, Bed, BigBed, Bam, VCF, and BCF files, which are distinguished thanks to their suffix (.wig, .bw, .bg, .bed, .bb, .bam, .vcf, .bcf respectively).");
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
puts("\tprogram = (iterator) | do (iterator) | (extraction) | (statistic)");
puts("\titerator = (in_filename) | (unary_operator) (iterator) | (binary_operator) (iterator) (iterator) | (reducer) (multiplex) | (setComparison) (multiplex_list) | print (output) (statistic)");
puts("\tunary_operator = unit | write (output) | write_bg (ouput) | smooth (int) | exp | ln | log (float) | pow (float) | offset (float) | scale (float) | gt (float) | lt (float) | default (float) | isZero | (statistic)");
puts("\toutput = (out_filename) | -");
puts("\tin_filename = *.wig | *.bw | *.bed | *.bb | *.bg | *.bam | *.vcf | *.bcf");
puts("\tstatistic = (statistic_function) (iterator)");
puts("\tstatistic_function = AUC | meanI | varI | minI | maxI | stddevI | CVI | pearson (iterator)");
puts("\tbinary_operator = diff | ratio | overlaps | apply (statistic) | fillIn");
puts("\treducer = cat | sum | product | mean | var | stddev | entropy | CV | median | min | max");
puts("\tsetComparison = ttest | ftest | wilcoxon");
puts("\tmultiplex_list = (multiplex) | (multiplex) : (multiplex_list)");
puts("\tmultiplex = (iterator_list) | map (unary_operator) (multiplex) | strict (multiplex)");
puts("\titerator_list = (iterator) | (iterator) : (iterator_list)");
puts("\textraction = profile (output) (int) (iterator) (iterator) | profiles (output) (int) (iterator) (iterator) | histogram (output) (width) (iterator_list) | mwrite (output) (multiplex) | mwrite_bg (output) (multiplex)");
puts("\t\t| apply_paste (out_filename) (statistic) (bed_file) (iterator)");

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

static WiggleIterator ** readIteratorList(int * count, bool * strict);

static WiggleIterator ** readMappedIteratorList(int * count, bool * strict) {
	char * token = needNextToken();
	WiggleIterator ** iters;
	int i;

	if (strcmp(token, "unit") == 0) {
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = UnitWiggleIterator(iters[i]);
	} else if (strcmp(token, "smooth") == 0) {
		int width = atoi(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = SmoothWiggleIterator(iters[i], width);
	} else if (strcmp(token, "exp") == 0) {
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = NaturalExpWiggleIterator(iters[i]);
	} else if (strcmp(token, "ln") == 0) {
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = NaturalLogWiggleIterator(iters[i]);
	} else if (strcmp(token, "log") == 0) {
		double base = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = LogWiggleIterator(iters[i], base);
	} else if (strcmp(token, "pow") == 0) {
		double base = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = PowerWiggleIterator(iters[i], base);
	} else if (strcmp(token, "scale") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = ScaleWiggleIterator(iters[i], scalar);
	} else if (strcmp(token, "offset") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = ShiftWiggleIterator(iters[i], scalar);
	} else if (strcmp(token, "gt") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = HighPassFilterWiggleIterator(iters[i], scalar);
	} else if (strcmp(token, "lt") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = HighPassFilterWiggleIterator(ScaleWiggleIterator(iters[i], -1), -scalar);
	} else if (strcmp(token, "default") == 0) {
		double scalar = atof(needNextToken());
		iters = readIteratorList(count, strict);
		for (i = 0; i < *count; i++)
			iters[i] = DefaultValueWiggleIterator(iters[i], scalar);
	} else {
		fprintf(stderr, "Unary function unkown: %s\n", token);
		exit(1);
	}

	return iters;
}

static WiggleIterator ** readIteratorListToken(int * count, bool * strict, char * token) {
	if (strcmp(token, "map") == 0)
		return readMappedIteratorList(count, strict);
	else if (strcmp(token, "strict") == 0) {
		*strict = true;
		return readIteratorList(count, strict);
	} else
		return readFileList(count, token);
}

static WiggleIterator ** readIteratorList(int * count, bool * strict) {
	return readIteratorListToken(count, strict, needNextToken());
}

static void noTokensLeft() {
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
}

static WiggleIterator ** readLastIteratorList(int * count) {
	bool strict = false;
	WiggleIterator ** res = readIteratorList(count, &strict);
	noTokensLeft();
	return res;
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

typedef WiggleIterator * (*statisticCreator)(WiggleIterator *);

static statisticCreator * readStatisticList(char ** token, int * count) {
	int maxLength = 8;
	*count = 0;
	statisticCreator * statistics = (statisticCreator *) calloc(maxLength, sizeof(statisticCreator));

	while (true) {
		if (strcmp(*token, "AUC") == 0)
			statistics[(*count)++] = &AUCIntegrator;
		else if (strcmp(*token, "meanI") == 0)
			statistics[(*count)++] = &MeanIntegrator;
		else if (strcmp(*token, "varI") == 0)
			statistics[(*count)++] = &VarianceIntegrator;
		else if (strcmp(*token, "stddevI") == 0)
			statistics[(*count)++] = &StandardDeviationIntegrator;
		else if (strcmp(*token, "CVI") == 0)
			statistics[(*count)++] = &CoefficientOfVariationIntegrator;
		else if (strcmp(*token, "maxI") == 0)
			statistics[(*count)++] = &MaxIntegrator;
		else if (strcmp(*token, "minI") == 0)
			statistics[(*count)++] = &MinIntegrator;
		else
			break;

		if (*count == maxLength) {
			maxLength *= 2;
			statistics = realloc(statistics, maxLength * sizeof(void *));
		}

		*token = needNextToken();
	}

	if (*count == 0) {
		fprintf(stderr, "Name of function to be applied unrecognized: %s\n", *token);
		exit(1);
	}

	return statistics;
}

static Multiplexer * readApply() {
	char * token = needNextToken();
	bool strict = true;
	int count;
	statisticCreator * statistics = readStatisticList(&token, &count);

	if (strcmp(token, "fillIn") == 0) {
		strict = false;
		token = needNextToken();
	}

	WiggleIterator * regions = readIteratorToken(token);
	WiggleIterator * data = readIterator();

	return ApplyMultiplexer(regions, statistics, count, data, strict);
}


static Multiplexer * readMultiplexer();

static Multiplexer * readMultiplexerToken(char * token) {
	if (strcmp(token, "mwrite") == 0) {
		FILE * file = readOutputFilename();
		return TeeMultiplexer(readMultiplexer(), file, false, holdFire);
	} else if (strcmp(token, "mwrite_bg") == 0) {
		FILE * file = readOutputFilename();
		return TeeMultiplexer(readMultiplexer(), file, true, holdFire);
	} else if (strcmp(token, "apply") == 0) {
		return readApply();
	} else {
		int count = 0; 
		bool strict = false;
		WiggleIterator ** iters = readIteratorListToken(&count, &strict, token);
		return newMultiplexer(iters, count, strict);
	}
}

static Multiplexer * readLastMultiplexerToken(char * token) {
	Multiplexer * res = readMultiplexerToken(token);
	noTokensLeft();
	return res;
}

static Multiplexer * readMultiplexer() {
	char * token = needNextToken();
	return readMultiplexerToken(token);
}

static Multiplexer ** readMultiplexerList(int * count) {
	size_t buffer_size = 8;
	char * token;
	Multiplexer ** multis = (Multiplexer **) calloc(buffer_size, sizeof(Multiplexer*));
	*count = 0;

	for (token = nextToken(0,0); token != NULL && strcmp(token, ":"); token = nextToken(0,0)) {
		if (*count == buffer_size) {
			buffer_size *= 2;
			multis = (Multiplexer **) realloc(multis, buffer_size * sizeof(Multiplexer*));
		}
		multis[(*count)++] = readMultiplexerToken(token);
	}
	return multis;
}

static Multiset * readMultiset() {
	int count = 0; 
	Multiplexer ** multis = readMultiplexerList(&count);
	return newMultiset(multis, count);
}

static WiggleIterator * readTee() {
	FILE * file = readOutputFilename();
	return TeeWiggleIterator(readIterator(), file, false, holdFire);
}

static WiggleIterator * readBGTee() {
	FILE * file = readOutputFilename();
	return TeeWiggleIterator(readIterator(), file, true, holdFire);
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

static WiggleIterator * readLt() {
	double cutoff = atof(needNextToken());
	return HighPassFilterWiggleIterator(ScaleWiggleIterator(readIterator(), -1), -cutoff);
}

static WiggleIterator * readDefault() {
	double value = atof(needNextToken());
	return DefaultValueWiggleIterator(readIterator(), value);
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

static WiggleIterator * readFillIn() {
	return FillInReduction(readMultiplexer());
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

static WiggleIterator * readPrint() {
	FILE * file = readOutputFilename();
	WiggleIterator * wi = readIterator();
	return PrintStatisticsWiggleIterator(wi, file);
}

static WiggleIterator * readDifference() {
	WiggleIterator ** iters = calloc(2, sizeof(WiggleIterator *));
	bool strict = false;
	char * token = needNextToken();
	if (strcmp(token, "strict")==0) {
		strict = true;
		token = needNextToken();
	}
	iters[0] = readIteratorToken(token);
	iters[1] = ScaleWiggleIterator(readIterator(), -1);
	return SumReduction(newMultiplexer(iters, 2, strict));
}

static WiggleIterator * readRatio() {
	WiggleIterator ** iters = calloc(2, sizeof(WiggleIterator *));
	bool strict = false;
	char * token = needNextToken();
	if (strcmp(token, "strict")==0) {
		strict = true;
		token = needNextToken();
	}
	iters[0] = readIteratorToken(token);
	iters[1] = PowerWiggleIterator(readIterator(), -1);
	return ProductReduction(newMultiplexer(iters, 2, strict));
}

static WiggleIterator * readSeek() {
	char * chrom = needNextToken();
	int start = atoi(needNextToken());
	int finish = atoi(needNextToken());
	holdFire = true;

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

static WiggleIterator * readFTest() {
	return FTestReduction(readMultiset());
}

static WiggleIterator * readMWUTest() {
	Multiplexer ** multis = calloc(2, sizeof(Multiplexer *));
	multis[0] = readMultiplexer();
	multis[1] = readMultiplexer();
	return MWUReduction(newMultiset(multis, 2));
}

static WiggleIterator * readMeanIntegrator() {
	return MeanIntegrator(readIterator());
}

static WiggleIterator * readMaxIntegrator() {
	return MaxIntegrator(readIterator());	
}

static WiggleIterator * readMinIntegrator() {
	return MinIntegrator(readIterator());
}

static WiggleIterator * readVarianceIntegrator() {
	return VarianceIntegrator(readIterator());	
}

static WiggleIterator * readStandardDeviationIntegrator() {
	return StandardDeviationIntegrator(readIterator());	
}

static WiggleIterator * readCoefficientOfVariationIntegrator() {
	return CoefficientOfVariationIntegrator(readIterator());	
}

static WiggleIterator * readPearson() {
	WiggleIterator * iter1 = readIterator();
	WiggleIterator * iter2 = readLastIterator();
	return PearsonIntegrator(iter1, iter2);
}

static WiggleIterator * readAUC() {
	return AUCIntegrator(readLastIterator());
}

static WiggleIterator * readIsZero() {
	return IsZero(readIterator());	
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
	if (strcmp(token, "print") == 0)
		return readPrint();
	if (strcmp(token, "sum") == 0)
		return readSum();
	if (strcmp(token, "fillIn") == 0)
		return readFillIn();
	if (strcmp(token, "mult") == 0)
		return readProduct();
	if (strcmp(token, "diff") == 0)
		return readDifference();
	if (strcmp(token, "ratio") == 0)
		return readRatio();
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
	if (strcmp(token, "lt") == 0)
		return readLt();
	if (strcmp(token, "default") == 0)
		return readDefault();
	if (strcmp(token, "overlaps") == 0)
		return readOverlap();
	if (strcmp(token, "ttest") == 0)
		return readTTest();
	if (strcmp(token, "ftest") == 0)
		return readFTest();
	if (strcmp(token, "wilcoxon") == 0)
		return readMWUTest();
	if (strcmp(token, "AUC") == 0)
		return readAUC();
	if (strcmp(token, "meanI") == 0)
		return readMeanIntegrator();
	if (strcmp(token, "maxI") == 0)
		return readMaxIntegrator();
	if (strcmp(token, "minI") == 0)
		return readMinIntegrator();
	if (strcmp(token, "varI") == 0)
		return readVarianceIntegrator();
	if (strcmp(token, "stddevI") == 0)
		return readStandardDeviationIntegrator();
	if (strcmp(token, "CVI") == 0)
		return readCoefficientOfVariationIntegrator();
	if (strcmp(token, "pearson") == 0) 
		return readPearson();
	if (strcmp(token, "isZero") == 0) 
		return readIsZero();
	if (strcmp(token, "apply") == 0) 
		return SelectReduction(readApply(), 0);

	return SmartReader(token, holdFire);

}

static void readProfile() {
	FILE * file = readOutputFilename();

	int width = atoi(needNextToken());
	WiggleIterator * regions = readIterator();
	WiggleIterator * wig = readLastIterator();
	Multiplexer * profiles = ProfileMultiplexer(regions, width, wig);
	double * profile = calloc(width, sizeof(double));

	for (; !profiles->done; popMultiplexer(profiles))
		addProfile(profile, profiles->values, width);

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
	Multiplexer * profiles;

	for (profiles = ProfileMultiplexer(regions, width, wig); !profiles->done; popMultiplexer(profiles)) {
		fprintf(file, "%s\t%i\t%i\t", profiles->chrom, profiles->start, profiles->finish);
		fprintfProfile(file, profiles->values, width);
	}

	fclose(file);
}

static void readHistogram() {
	FILE * file = readOutputFilename();
	int width = atoi(needNextToken());
	int count = 0; 
	WiggleIterator ** iters = readLastIteratorList(&count);
	Histogram * hist = histogram(iters, count, width);	
	print_histogram(hist, file);
	fclose(file);
}

static Multiplexer * readApplyPaste() {
	FILE * outfile = readOutputFilename();
	bool strict = true;
	int count;
	char * token = needNextToken();
	statisticCreator * statistics = readStatisticList(&token, &count);

	if (strcmp(token, "fillIn") == 0) {
		strict = false;
		token = needNextToken();
	}

	char * infilename = token;
	FILE * infile = fopen(token, "r");
	if (!infile) {
	       fprintf(stderr, "Could not open %s.\n", token);
	       exit(1);
	}

	return PasteMultiplexer(ApplyMultiplexer(SmartReader(infilename, holdFire), statistics, count, readLastIterator(), strict), infile, outfile, false);
}

void rollYourOwn(int argc, char ** argv) {
	char * token = nextToken(argc, argv);
	if (strcmp(token, "do") == 0)
		runWiggleIterator(readLastIterator());
	else if (strncmp(token, "write", 5) == 0)
		runWiggleIterator(readLastIteratorToken(token));
	else if (strncmp(token, "mwrite", 6) == 0)
		runMultiplexer(readLastMultiplexerToken(token));
	else if (strcmp(token, "apply_paste") == 0)
		runMultiplexer(readApplyPaste());
	else if (strcmp(token, "histogram") == 0)
		readHistogram();
	else if (strcmp(token, "profile") == 0)
		readProfile();
	else if (strcmp(token, "profiles") == 0)
		readProfiles();
	else if (strcmp(token, "print") == 0)
		runWiggleIterator(readLastIteratorToken(token));
	else if (strcmp(token, "AUC") == 0 || strcmp(token, "meanI") == 0 || strcmp(token, "varI") == 0 || strcmp(token, "stddevI") == 0 || strcmp(token, "CVI") == 0 || strcmp(token, "maxI") == 0 || strcmp(token, "minI") == 0 || strcmp(token, "pearson") == 0)
		runWiggleIterator(PrintStatisticsWiggleIterator(readLastIteratorToken(token), stdout));
	else
		toStdout(readLastIteratorToken(token), false, false);	
}
