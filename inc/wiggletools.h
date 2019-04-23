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

#ifndef _WIGGLETOOLS_DEF_
#define _WIGGLETOOLS_DEF_

#ifndef bool
#define bool char
#define true 1
#define false 0
#endif

#include <stdio.h>

typedef struct wiggleIterator_st WiggleIterator;
typedef struct multiplexer_st Multiplexer;
typedef struct multiset_st Multiset;
typedef struct histogram_st Histogram;

// Creators
WiggleIterator * SmartReader (char *, bool);
WiggleIterator * CatWiggleIterator (char **, int);
// Secondary creators (to force file format recognition if necessary)
WiggleIterator * WiggleReader (char *);
WiggleIterator * BigWiggleReader (char *, bool);
WiggleIterator * BedReader (char *);
WiggleIterator * BigBedReader (char *, bool);
WiggleIterator * BamReader (char *, bool);
WiggleIterator * SamReader (char *);
WiggleIterator * VcfReader (char *);
WiggleIterator * BcfReader (char *, bool);

// Generic class functions 
void seek(WiggleIterator *, const char *, int, int);

// Algebraic operations on iterators
	
	// Unary
WiggleIterator * UnitWiggleIterator (WiggleIterator *);
WiggleIterator * CoverageWiggleIterator (WiggleIterator *);
WiggleIterator * UnionWiggleIterator (WiggleIterator *);
WiggleIterator * NonOverlappingWiggleIterator (WiggleIterator *);
WiggleIterator * AbsWiggleIterator (WiggleIterator * );
WiggleIterator * NaturalLogWiggleIterator (WiggleIterator *);
WiggleIterator * NaturalExpWiggleIterator (WiggleIterator *);
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator * );
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator *);
WiggleIterator * OverlapWiggleIterator(WiggleIterator *, WiggleIterator *);
WiggleIterator * TrimWiggleIterator(WiggleIterator *, WiggleIterator *);
WiggleIterator * NoverlapWiggleIterator(WiggleIterator *, WiggleIterator *);
WiggleIterator * NearestWiggleIterator(WiggleIterator *, WiggleIterator *);
WiggleIterator * IsZero(WiggleIterator *);
WiggleIterator * Floor(WiggleIterator *);
	// Scalar operations
WiggleIterator * ScaleWiggleIterator (WiggleIterator *, double);
WiggleIterator * ShiftWiggleIterator(WiggleIterator *, double);
WiggleIterator * PowerWiggleIterator (WiggleIterator *, double);
WiggleIterator * LogWiggleIterator (WiggleIterator * , double);
WiggleIterator * ExpWiggleIterator (WiggleIterator *, double);
WiggleIterator * DefaultValueWiggleIterator(WiggleIterator *, double);
WiggleIterator * HighPassFilterWiggleIterator(WiggleIterator *, double);
WiggleIterator * SmoothWiggleIterator(WiggleIterator * i, int);
WiggleIterator * ExtendWiggleIterator(WiggleIterator * i, int);

// Sets of iterators 
Multiplexer * newMultiplexer(WiggleIterator **, int, bool);

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
WiggleIterator * FillInReduction( Multiplexer * );

// Sets of sets iterators 
Multiset * newMultiset(Multiplexer **, int);

// Reduction operators on sets of sets:
WiggleIterator * TTestReduction(Multiset *);
WiggleIterator * FTestReduction(Multiset *);
WiggleIterator * MWUReduction(Multiset *);

// Output
void toFile (WiggleIterator *, char *, bool, bool);
void toStdout (WiggleIterator *, bool, bool);
WiggleIterator * TeeWiggleIterator(WiggleIterator *, FILE *, bool, bool);
void runWiggleIterator(WiggleIterator * );
Multiplexer * TeeMultiplexer(Multiplexer *, FILE *, bool, bool);
void toStdoutMultiplexer (Multiplexer *, bool, bool);
void runMultiplexer(Multiplexer * );
WiggleIterator * PrintStatisticsWiggleIterator(WiggleIterator * i, FILE * file);

// Statistics
// 	Unary
WiggleIterator * AUCIntegrator (WiggleIterator *);
WiggleIterator * MeanIntegrator (WiggleIterator *);
WiggleIterator * MinIntegrator (WiggleIterator *);
WiggleIterator * MaxIntegrator (WiggleIterator *);
WiggleIterator * VarianceIntegrator (WiggleIterator *);
WiggleIterator * StandardDeviationIntegrator (WiggleIterator *);
WiggleIterator * CoefficientOfVariationIntegrator (WiggleIterator *);
WiggleIterator * NDPearsonIntegrator(Multiset *);
void regionProfile(WiggleIterator *, double *, int, int, bool);
void addProfile(double *, double *, int);
//	Binary 
WiggleIterator * PearsonIntegrator (WiggleIterator * , WiggleIterator * );
//	Histograms
Histogram * histogram(WiggleIterator **, int, int);
void normalize_histogram(Histogram *);
void print_histogram(Histogram *, FILE *);

// Regional statistics
Multiplexer * ApplyMultiplexer(WiggleIterator *, WiggleIterator * (**statistics)(WiggleIterator *), int count, WiggleIterator *, bool strict);
Multiplexer * ProfileMultiplexer(WiggleIterator *, int, WiggleIterator *);
Multiplexer * PasteMultiplexer(Multiplexer *,  FILE *, FILE *, bool);

// Cleaning up
void destroyWiggleIterator (WiggleIterator *);

// Big file params
void setMaxBlocks(int);
void setMaxHeadStart(int);

// Command line parser
void rollYourOwn(int argc, char ** argv);
void printHelp();

#endif
