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

// Local header
#include "bigFileReader.h"

int MAX_BLOCKS = 100;
int MAX_HEAD_START = 300;

static BlockData * createBlockData(char * chrom, struct fileOffsetSize * block, char * blockBuf, size_t uncompressBufSize) {
	BlockData * new = (BlockData * ) calloc(1, sizeof(BlockData));
	new->chrom = chrom;
	if (uncompressBufSize > 0) {
		new->duplicate = true;
		new->uncompressBuf = (char *) needLargeMem(uncompressBufSize);
		int uncSize = zUncompress(blockBuf, block->size, new->uncompressBuf, uncompressBufSize);
		new->blockEnd = new->uncompressBuf + uncSize;
	} else {
		new->duplicate = false;
		new->uncompressBuf = blockBuf;
		new->blockEnd = blockBuf + block->size;
	}
	pthread_mutex_init(&new->proceed, NULL);
	pthread_mutex_lock(&new->proceed);
	return new;
}

void destroyBlockData(BlockData * data) {
	if (data->duplicate)
		free(data->uncompressBuf);
	pthread_mutex_destroy(&data->proceed);
	free(data);
}

void enterBlock(BigFileReaderData * data) {
	data->blockPt = data->blockData->uncompressBuf;
	data->chrom = data->blockData->chrom;
}

void goToNextBlock(BigFileReaderData * data) {
	BlockData * prevBlockData = data->blockData;

	// Check whether allowed to step forward
	pthread_mutex_lock(&(data->blockData->proceed));
	data->blockData = data->blockData->next;
	pthread_mutex_destroy(&prevBlockData->proceed);
	destroyBlockData(prevBlockData);

	// Signal freed memory
	pthread_mutex_lock(&data->count_mutex);
	data->blockRuns--;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
}

static void waitForPermission(BigFileReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	if (data->blockRuns > MAX_HEAD_START)
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	data->blockRuns++;
	pthread_mutex_unlock(&data->count_mutex);
}

static void freeLastBarrier(BigFileReaderData * data) {
	if (data->lastBlockData) {
		pthread_mutex_unlock(&(data->lastBlockData->proceed));
	} else { 
		pthread_mutex_lock(&data->proceed_mutex);
		pthread_cond_signal(&data->proceed_cond);
		pthread_mutex_unlock(&data->proceed_mutex);
	}
}

static void * downloadBlockRun(BigFileReaderData * data, char * chrom, struct fileOffsetSize * firstBlock, struct fileOffsetSize * afterBlock, bits64 mergedSize) {
	char * mergedBuf, *blockBuf;
	struct fileOffsetSize * block;

	udcSeek(data->udc, firstBlock->offset);
	blockBuf = mergedBuf = (char *) needLargeMem(mergedSize);
	udcMustRead(data->udc, mergedBuf, mergedSize);

	for (block = firstBlock; block != afterBlock; block = block->next) {
		waitForPermission(data);
		BlockData * new = createBlockData(chrom, block, blockBuf, data->bwf->uncompressBufSize);
		if (data->lastBlockData) {
			data->lastBlockData->next = new;
		} else {
			data->blockData = new;
		}

		freeLastBarrier(data);

		// Move on
		data->lastBlockData = new;
		blockBuf += block->size;
	}

	freeMem(mergedBuf);
	return NULL;
}

static void downloadBigRegion(BigFileReaderData * data, char * chrom, int start, int finish) {
	struct fileOffsetSize *blockList, *block, *beforeGap, *afterGap;
	int blockCounter;
	bits64 mergedSize;

	blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, chrom, start, finish, NULL);

	for (block = blockList; block; block=block->next) {
		/* Read contiguous blocks into mergedBuf. */
		fileOffsetSizeFindGap(block, &beforeGap, &afterGap);

		// Little hack to limit the number of blocks read at any time
		struct fileOffsetSize * blockPtr, * prevBlock;
		blockCounter = 0;
		prevBlock = block;
		for (blockPtr = block; blockPtr != afterGap && blockCounter < MAX_BLOCKS; blockPtr = blockPtr->next) {
			blockCounter++;
			prevBlock = blockPtr;
		}
		if (blockCounter == MAX_BLOCKS) {
			beforeGap = prevBlock;
			afterGap = blockPtr;
		}

		mergedSize = beforeGap->offset + beforeGap->size - block->offset;

		downloadBlockRun(data, chrom, block, afterGap, mergedSize);
	}

	slFreeList(blockList);
}

static void downloadFullGenome(BigFileReaderData * data) {
	struct bbiChromInfo *chromList, *chrom;
	chromList = bbiChromList(data->bwf);

	for (chrom = chromList; chrom; chrom = chrom->next) 
		downloadBigRegion(data, chrom->name, 0, chrom->size);
	// TODO free chromList memory... yes but labels lost! Need to be copied out first....

}

static void * downloadBigFile(void * args) {
	BigFileReaderData * data = (BigFileReaderData *) args;

	if (!data->chrom)
		downloadFullGenome(data);
	else 
		downloadBigRegion(data, data->chrom, data->start, data->stop);

	freeLastBarrier(data);
	return NULL;
}

void BigFileReaderCloseFile(BigFileReaderData * data) {
	//bbiChromInfoFreeList(&(data->chromList));
	bbiFileClose(&(data->bwf));
}

void launchDownloader(BigFileReaderData * data) {
	pthread_mutex_init(&data->proceed_mutex, NULL);
	pthread_cond_init(&data->proceed_cond, NULL);
	pthread_mutex_init(&data->count_mutex, NULL);
	pthread_cond_init(&data->count_cond, NULL);

	// Launch download thread
	int err = pthread_create(&data->downloaderThreadID, NULL, &downloadBigFile, data);
	if (err) {
		printf("Could not create new thread %i\n", err);
		abort();
	}
	pthread_detach(data->downloaderThreadID);

	// Wait for the first block to be available
	pthread_mutex_lock(&data->proceed_mutex);
	if (!data->blockData)
		pthread_cond_wait(&data->proceed_cond, &data->proceed_mutex);
	pthread_mutex_unlock(&data->proceed_mutex);
}

void killDownloader(BigFileReaderData * data) {
	if (data->downloaderThreadID)
		pthread_cancel(data->downloaderThreadID);

	pthread_mutex_destroy(&data->proceed_mutex);
	pthread_cond_destroy(&data->proceed_cond);
	pthread_mutex_destroy(&data->count_mutex);
	pthread_cond_destroy(&data->count_cond);

	while (data->blockData) {
		BlockData * prevData = data->blockData;
		data->blockData = data->blockData->next;
		destroyBlockData(prevData);
	}

	data->lastBlockData = NULL;
	data->blockRuns = 0;
}
