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

#ifndef _BIG_WIGGLE_READER_H_
#define _BIG_WIGGLE_READER_H_

#include <stdio.h>
#include <pthread.h>

#include "wiggleIterator.h"

// Kent library headers
#include "common.h"
#include "bigBed.h"
#include "bbiFile.h"
#include "zlibFace.h"
#include "bbiFile.h"
#include "bigWig.h"
#include "udc.h"
#include "bwgInternal.h"

typedef struct blockData_st {
	char *uncompressBuf;
	char *blockEnd;
	struct blockData_st * next;
	boolean duplicate;
	char *chrom;
} BlockData;

typedef struct bigFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int start, stop;
	pthread_t downloaderThreadID;

	// Output of downloader
	struct bbiFile* bwf;
	struct udcFile *udc;
	BlockData * blockData, *lastBlockData;
	pthread_mutex_t count_mutex;
	pthread_cond_t count_cond;
	int blockCount;

	// Items within block
	struct bwgSectionHead head;
	boolean isSwapped;
	char *blockPt;
	bits16 i;

} BigFileReaderData;

void launchDownloader(BigFileReaderData * data);
void killDownloader(BigFileReaderData * data);
void destroyBlockData(BlockData * data);
void enterBlock(BigFileReaderData * data);
void goToNextBlock(BigFileReaderData * data);
#endif
