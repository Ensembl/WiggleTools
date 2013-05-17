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
	puts("Help!");
	puts("Inputs:");
	puts("\tThe program takes in Wig and BigWig files, which are distinguished thanks to their suffix (.wig or .bw respectively). There are several modes to run wiggletools:");
	puts("");
	puts("Outputs:");
	puts("\tThe program outputs a flat file wiggle file in stdout.");
	puts("");
	puts("\tParameters:");
	puts("\t// Unary operators");
	puts("\twiggletools exp file");
	puts("\twiggletools log file");
	puts("\t");
	puts("\t// Operators between a signal a scalar");
	puts("\twiggletools scale file factor");
	puts("\twiggletools pow file exponent");
	puts("\twiggletools exp file radix");
	puts("\twiggletools log file base");
	puts("\t");
	puts("\t// Binary operators between two signal files");
	puts("\twiggletools add file1 file2");
	puts("\twiggletools mult file1 file2");
}

int main(int argc, char ** argv) {
	if (argc < 2 || strcmp(argv[1], "help") == 0) {
		printHelp();
		return 0;
	} else if (strcmp(argv[1], "add") == 0) {
		toStdout(SumWiggleIterator(WigOrBigWigReader(argv[2]), WigOrBigWigReader(argv[3])));
	} else if (strcmp(argv[1], "scale") == 0) {
		toStdout(ScaleWiggleIterator(WigOrBigWigReader(argv[2]), atoi(argv[3])));
	} else if (strcmp(argv[1], "mult") == 0) {
		toStdout(ProductWiggleIterator(WigOrBigWigReader(argv[2]), WigOrBigWigReader(argv[3])));
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

	return 1;
}

