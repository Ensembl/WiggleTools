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
	puts("wiggletools [opts] 'command'");
	puts("");
	puts("Options:");
	puts("\t-maxBlocks n\t: set the max number of blocks read at once in a BigFile reader");
	puts("\t-maxHeadStart n\t: set the max number of decompressed blocks in a BigFile reader");
	puts("\t--help");
	puts("");
	puts("Command grammar:");
	puts("\tcommand:\t\titerator|statistic iterator|apply|comparison");
	puts("\titerator:\t\tfilename|unary_operation|binary_operation|reduction|output_operation");
	puts("\tunary_operation:\tunary_operator iterator");
	puts("\tunary_operator:\t\t'unit'");
	puts("\tbinary_operation:\tbinary_operator iterator");
	puts("\tbinary_operator:\t'diff'");
	puts("\treduction:\t\treduction_operator multiplexer");
	puts("\treduction_operator:\t'cat'|'add'|'product'|'mean'|'variance'|'stddev'|'median'|'min'|'max'");
	puts("\tmultiplexer:\t\tfilename_list|map");
	puts("\tfilename list:\t\tfilename1 filename2 ... ; ");
	puts("\tmap:\t\t\t'map' unary_operator multiplexer");
	puts("\toutput_operation:\toutput_operator filename iterator");
	puts("\toutput_operator:\t'write'|'writeb'|'print'");
	puts("\tstatistic:\t\t'AUC'|'mean'");
	puts("\tapply:\t\t\t'apply' statistic iterator iterator");
	puts("\tcomparison:\t\tcomparator iterator iterator");
	puts("\tcomparator:\t\t'pearson'");
}

int main(int argc, char ** argv) {
	int i=0;
	if (argc < 2 || strcmp(argv[1], "--help") == 0) {
		printHelp();
		return 0;
	}

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-maxBlocks") == 0) {
			i++;
			int value;
			sscanf(argv[i], "%i", &value);
			setMaxBlocks(value);
		} else if (strcmp(argv[i], "-maxHeadStart") == 0) {
			i++;
			int value;
			sscanf(argv[i], "%i", &value);
			setMaxHeadStart(value);
		} else {
			rollYourOwn(argv[i]);
		}
	}

	return 0;
}

