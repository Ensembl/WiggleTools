// Copyright [1999-2016] EMBL-European Bioinformatics Institute
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

#ifndef _BIG_WIGGLE_READER_H_
#define _BIG_WIGGLE_READER_H_

#include <stdio.h>

#include "wiggleIterator.h"
#include "bufferedReader.h"

// Kent library headers
#include "common.h"
#include "bigBed.h"
#include "bbiFile.h"
#include "zlibFace.h"
#include "bbiFile.h"
#include "bigWig.h"
#include "udc.h"
#include "bwgInternal.h"

typedef struct bigFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int start, stop;
	bool (*readBuffer)(struct bigFileReaderData_st *);

	// Output of downloader
	BufferedReaderData * bufferedReaderData;

	// BigFile variables
	struct bbiFile* bwf;
	struct udcFile *udc;
	boolean isSwapped;

	// Buffer data
	char *uncompressBuf;
	char *blockEnd;
} BigFileReaderData;

void openBigFile(BigFileReaderData * data);
void * downloadBigFile(void * data);
void BigFileReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish);
void BigFileReaderPop(WiggleIterator * wi);
void killDownloader(BigFileReaderData * data);
#endif
