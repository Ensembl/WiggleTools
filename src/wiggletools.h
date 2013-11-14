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

// Generic class functions 
void seek(WiggleIterator *, const char *, int, int);

// Algebraic operations on iterators
	
	// Unary
WiggleIterator * UnitWiggleIterator (WiggleIterator *);
WiggleIterator * UnionWiggleIterator (WiggleIterator *);
WiggleIterator * AbsWiggleIterator (WiggleIterator * );
WiggleIterator * NaturalLogWiggleIterator (WiggleIterator *);
WiggleIterator * NaturalExpWiggleIterator (WiggleIterator *);
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator * );
WiggleIterator * HighPassFilterWiggleIterator(WiggleIterator *, double);
WiggleIterator * TestNonOverlappingWiggleIterator(WiggleIterator *);
	// Scalar operations
WiggleIterator * ScaleWiggleIterator (WiggleIterator *, double);
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
WiggleIterator * MedianReduction ( Multiplexer * );

// Sets of sets iterators 
Multiplexer * newMultiset(Multiplexer **, int);

// Reduction operators on sets of sets:
WiggleIterator * TTestReduction(Multiset *);
WiggleIterator * MWUReduction(Multiset *);

// Output
void toFile (WiggleIterator *, char *);
void toBinaryFile (WiggleIterator *, char *);
void toStdout (WiggleIterator *);
WiggleIterator * BinaryTeeWiggleIterator(WiggleIterator *, FILE *);
WiggleIterator * TeeWiggleIterator(WiggleIterator *, FILE *);
void runWiggleIterator(WiggleIterator * );

WiggleIterator * ApplyWiggleIterator(WiggleIterator * regions, double (*statistic)(WiggleIterator *), WiggleIterator * dataset);
WiggleIterator * PasteWiggleIterator(WiggleIterator * i, FILE * infile, FILE * outfile);

// Statistics
// 	Unary
double AUC (WiggleIterator *);
double mean (WiggleIterator *);
double variance (WiggleIterator *);
double isZero(WiggleIterator * wi);
//	Binary 
double pearsonCorrelation (WiggleIterator * , WiggleIterator * );

// Regional statistics
WiggleIterator * apply(WiggleIterator * , double (*statistic)(WiggleIterator *), WiggleIterator *);

// Cleaning up
void destroyWiggleIterator (WiggleIterator *);

// Big file params
void setMaxBlocks(int);
void setMaxHeadStart(int);

// Command line parser
void rollYourOwn(char *);
void printHelp();

// Deprecated
Multiplexer * newStreamingMultiplexer(FILE * input);
Multiplexer * newIteratorMultiplexer(WiggleIterator *, int, int);
void streamWiggleIteratorAtIndex(FILE * , WiggleIterator * , int , int );
void streamMultiplexer(FILE *, Multiplexer *);

#endif
