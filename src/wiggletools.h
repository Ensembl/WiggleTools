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

#ifndef _WIGGLETOOLS_DEF_
#define _WIGGLETOOLS_DEF_

#ifndef true
typedef int bool;
#define true 1
#define false 0
#endif

#include <stdio.h>

typedef struct wiggleIterator_st WiggleIterator;
typedef struct multiplexer_st Multiplexer;
typedef struct multiset_st Multiset;
typedef struct histogram_st Histogram;

// Creators
WiggleIterator * SmartReader (char *);
WiggleIterator * CatWiggleIterator (char **, int);
// Secondary creators (to force file format recognition if necessary)
WiggleIterator * WiggleReader (char *);
WiggleIterator * BigWiggleReader (char *);
WiggleIterator * BedReader (char *);
WiggleIterator * BigBedReader (char *);
WiggleIterator * BamReader (char *);
WiggleIterator * BinaryFileReader (char *);
WiggleIterator * VcfReader (char *);
WiggleIterator * BcfReader (char *);

// Generic class functions 
void seek(WiggleIterator *, const char *, int, int);

// Algebraic operations on iterators
	
	// Unary
WiggleIterator * UnitWiggleIterator (WiggleIterator *);
WiggleIterator * UnionWiggleIterator (WiggleIterator *);
WiggleIterator * NonOverlappingWiggleIterator (WiggleIterator *);
WiggleIterator * AbsWiggleIterator (WiggleIterator * );
WiggleIterator * NaturalLogWiggleIterator (WiggleIterator *);
WiggleIterator * NaturalExpWiggleIterator (WiggleIterator *);
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator * );
WiggleIterator * HighPassFilterWiggleIterator(WiggleIterator *, double);
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator *);
WiggleIterator * OverlapWiggleIterator(WiggleIterator * source, WiggleIterator * mask);
	// Scalar operations
WiggleIterator * ScaleWiggleIterator (WiggleIterator *, double);
WiggleIterator * ShiftWiggleIterator(WiggleIterator *, double);
WiggleIterator * PowerWiggleIterator (WiggleIterator *, double);
WiggleIterator * LogWiggleIterator (WiggleIterator * , double);
WiggleIterator * ExpWiggleIterator (WiggleIterator *, double);
WiggleIterator * SmoothWiggleIterator(WiggleIterator * i, int width);

// Sets of iterators 
Multiplexer * newMultiplexer(WiggleIterator **, int);

// Reduction operators on sets

WiggleIterator * SelectReduction(Multiplexer *, int);
WiggleIterator * MaxReduction ( Multiplexer * );
WiggleIterator * MinReduction ( Multiplexer * );
WiggleIterator * SumReduction ( Multiplexer * );
WiggleIterator * ProductReduction ( Multiplexer * );
WiggleIterator * MeanReduction ( Multiplexer * );
WiggleIterator * VarianceReduction ( Multiplexer * );
WiggleIterator * StdDevReduction ( Multiplexer * );
WiggleIterator * EntropyReduction ( Multiplexer * );
WiggleIterator * CVReduction ( Multiplexer * );
WiggleIterator * MedianReduction ( Multiplexer * );

// Sets of sets iterators 
Multiset * newMultiset(Multiplexer **, int);

// Reduction operators on sets of sets:
WiggleIterator * TTestReduction(Multiset *);
WiggleIterator * MWUReduction(Multiset *);

// Output
void toFile (WiggleIterator *, char *, bool);
void toBinaryFile (WiggleIterator *, char *, bool);
void toStdout (WiggleIterator *, bool);
WiggleIterator * BinaryTeeWiggleIterator(WiggleIterator *, FILE *, bool);
WiggleIterator * TeeWiggleIterator(WiggleIterator *, FILE *, bool);
void runWiggleIterator(WiggleIterator * );

WiggleIterator * ApplyWiggleIterator(WiggleIterator * regions, double (*statistic)(WiggleIterator *), WiggleIterator * dataset);
WiggleIterator * ProfileWiggleIterator(WiggleIterator * regions, int width, WiggleIterator * dataset);
WiggleIterator * PasteWiggleIterator(WiggleIterator * i, FILE * infile, FILE * outfile);

// Statistics
// 	Unary
double AUC (WiggleIterator *);
double mean (WiggleIterator *);
double variance (WiggleIterator *);
double isZero(WiggleIterator * wi);
void regionProfile(WiggleIterator *, double *, int, int, bool);
void addProfile(double *, double *, int);
//	Binary 
double pearsonCorrelation (WiggleIterator * , WiggleIterator * );
//	Histograms
Histogram * histogram(WiggleIterator *, int);
void normalize_histogram(Histogram *);
void print_histogram(Histogram *, FILE *);

// Regional statistics
WiggleIterator * apply(WiggleIterator * , double (*statistic)(WiggleIterator *), WiggleIterator *);

// Cleaning up
void destroyWiggleIterator (WiggleIterator *);

// Big file params
void setMaxBlocks(int);
void setMaxHeadStart(int);

// Command line parser
void rollYourOwn(int argc, char ** argv);
void printHelp();

// Deprecated
Multiplexer * newStreamingMultiplexer(FILE * input);
Multiplexer * newIteratorMultiplexer(WiggleIterator *, int, int);
void streamWiggleIteratorAtIndex(FILE * , WiggleIterator * , int , int );
void streamMultiplexer(FILE *, Multiplexer *);

#endif
