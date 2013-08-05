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
