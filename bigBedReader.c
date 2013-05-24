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
#include "bigWiggleReader.h"

//////////////////////////////////////////////////////
// Compressed File Reader
//////////////////////////////////////////////////////

static void * unzipBlockData(void * args) {
	BlockData * data = (BlockData *) args;
	data->uncompressBuf = (char *) needLargeMem(data->bwf->uncompressBufSize);
	int uncSize = zUncompress(data->blockBuf, data->block->size, data->uncompressBuf, data->bwf->uncompressBufSize);
	data->blockEnd = data->uncompressBuf + uncSize;
	data->done = true;

	// Communicate that the slot is free
	pthread_mutex_lock(data->count_mutex);
	(*(data->count))++;
	pthread_cond_signal(data->count_threshold_cv);
	pthread_mutex_unlock(data->count_mutex);

	return NULL;
}

static BlockData * createBlockData(BigWiggleReaderData * data, struct fileOffsetSize * block, pthread_mutex_t * count_mutex, pthread_cond_t * count_threshold_cv, int * count) {
	BlockData * new = (BlockData * ) calloc(1, sizeof(BlockData));
	new->blockBuf = data->blockBuf;
	new->block = block;
	new->bwf = data->bwf;
	new->blockEnd = NULL;
	new->done = false;
	new->count_mutex = count_mutex;
	new->count_threshold_cv = count_threshold_cv;
	new->count = count;
	pthread_mutex_init(&(new->proceed), NULL);

	// Launch thread...
	pthread_create(&(new->threadID), NULL, &unzipBlockData, new);
	return new;
}

static int findFreeThreadSlot(BlockData ** blockDatas) {
	int i;
	for (i = 0; i < THREADS; i++)
		if (!blockDatas[i] || blockDatas[i]->done)
			return i;

	exit(1);
}

static void waitForSlot(pthread_mutex_t * count_mutex, pthread_cond_t * count_threshold_cv, int * count) {
	pthread_mutex_lock(count_mutex);
	if (*count <= 0)
		pthread_cond_wait(count_threshold_cv, count_mutex);
	(*count)--;
	pthread_mutex_unlock(count_mutex);
}

static int getNextSlot(pthread_mutex_t * count_mutex, pthread_cond_t * count_threshold_cv, int * count, BlockData ** blockDatas) {
	waitForSlot(count_mutex, count_threshold_cv, count);
	return findFreeThreadSlot(blockDatas);
}

static void * createBlockDataList(void * args) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) args;
	struct fileOffsetSize * block;
	BlockData * last = NULL;
	BlockData * blockDatas[THREADS];
	int i;
	int count = THREADS;
	pthread_mutex_t count_mutex;
	pthread_cond_t count_threshold_cv;

	// Initialize mutexes
	pthread_mutex_init(&count_mutex, NULL);
	pthread_cond_init (&count_threshold_cv, NULL);
	
	// Initialize pointers to block data structs
	for (i = 0; i < THREADS; i++)
		blockDatas[i] = NULL;

	// Go through run of blocks
	for (block = data->block; block != data->afterGap; block = block->next) {
		i = getNextSlot(&count_mutex, &count_threshold_cv, &count, blockDatas);
		BlockData * new = blockDatas[i] = createBlockData(data, block, &count_mutex, &count_threshold_cv, &count);

		// Signal to reading thread that it is allowed to follow in list:
		pthread_mutex_lock(&(new->proceed));
		if (last) {
			last->next = new;
			pthread_mutex_unlock(&(last->proceed));
		} else { 
			data->blockData = new;
			pthread_barrier_wait(&(data->proceed));
		}

		// Move on
		last = new;
	}

	return NULL;
}

static void BigBedReaderEnterBlock(BigWiggleReaderData * data) {
	if (data->uncompress)
	    {
	    // Make sure data was unzipped
	    pthread_join(data->blockData->threadID, NULL);
	    data->blockPt = data->blockData->uncompressBuf;
	    data->blockEnd = data->blockData->blockEnd;
	    }
	else
	    {
	    data->blockPt = data->blockBuf;
	    data->blockEnd = data->blockPt + data->block->size;
	    }

	/* Deal with insides of block. */
	data->i = 0;
}

static void BigBedReaderEnterRunOfBlocks(BigWiggleReaderData * data) {
	/* Read contiguous blocks into mergedBuf. */
	struct fileOffsetSize *beforeGap;
	bits64 mergedSize;
	fileOffsetSizeFindGap(data->block, &beforeGap, &(data->afterGap));
	mergedSize = beforeGap->offset + beforeGap->size - data->block->offset;
	udcSeek(data->udc, data->block->offset);
	data->mergedBuf = (char *) needLargeMem(mergedSize);
	udcMustRead(data->udc, data->mergedBuf, mergedSize);
	data->blockBuf = data->mergedBuf;
	// If need to unzip, launch threads ahead of reader
	if (data->uncompress) {
		pthread_t threadID;

		pthread_barrier_init(&(data->proceed), NULL, 2);
		pthread_create(&threadID, NULL, &createBlockDataList, data);
		pthread_barrier_wait(&(data->proceed));
		pthread_barrier_destroy(&(data->proceed));
	}
	BigBedReaderEnterBlock(data); 
}

void BigBedReaderEnterChromosome(BigWiggleReaderData * data) {
	data->block = data->blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, data->chrom->name, 0, data->chrom->size, NULL);
	BigBedReaderEnterRunOfBlocks(data);
}

static void BigBedReaderCloseFile(BigWiggleReaderData * data) {
	// Because strings are passed by reference instead of pointers deallocating
	// this linked list would mess up objects for iterators downstream
	//bbiChromInfoFreeList(&(data->chromList));
	bbiFileClose(&(data->bwf));
}

static void BigBedReaderGoToNextChromosome(BigWiggleReaderData * data) {
	slFreeList(&(data->blockList));
	if ((data->chrom = data->chrom->next))
		BigBedReaderEnterChromosome(data);
	else 
		BigBedReaderCloseFile(data);
}

static void BigBedReaderGoToNextRunOfBlocks(BigWiggleReaderData * data) {
	freeMem(data->mergedBuf);
	if (data->block)
		BigBedReaderEnterRunOfBlocks(data);
	else if (data->chrom) 
		BigBedReaderGoToNextChromosome(data);
}

void BigBedReaderGoToNextBlock(BigWiggleReaderData * data) {
	data->blockBuf += data->block->size;
	if ((data->block = data->block->next) == data->afterGap)
		BigBedReaderGoToNextRunOfBlocks(data);
	else {
		if (data->uncompress) {
			BlockData * blockData = data->blockData;
			pthread_mutex_lock(&(blockData->proceed));
			data->blockData = blockData->next;
			freeMem(blockData->uncompressBuf);
			pthread_mutex_destroy(&(blockData->proceed));
			free(blockData);
		}
		BigBedReaderEnterBlock(data);
	}
}

static void BigBedReaderPop2(WiggleIterator * wi) {
	BigWiggleReaderData * data;

	if (wi->done)
		return;

	data = (BigWiggleReaderData*) wi->data;

	if (!data->chrom) {
		// Passive agressive indicator that the iterator should be closed
		// (avoids passing the WiggleIterator reference needlessly across 
		// all functions).
		wi->done = true;
		return;
	} else {
		wi->chrom = data->chrom->name;
	}

	/* Read next record into local variables. */
	memReadBits32(&data->blockPt, data->isSwapped);	// Read and discard chromId
	wi->start = memReadBits32(&data->blockPt, data->isSwapped);
	wi->finish = memReadBits32(&data->blockPt, data->isSwapped) + 1;
	// Skip boring stuff...
	for (;;)
		if (*(data->blockPt++) <= 0)
			break;

        if (data->blockPt == data->blockEnd)
	    BigBedReaderGoToNextBlock(data);
}

void BigBedReaderPop(WiggleIterator * wi) {
	int finish = wi->finish;
	char chrom[1000];

	strcpy(chrom, wi->chrom);

	while (wi->chrom && !wi->done) {
		BigBedReaderPop2(wi);
		if (chrom[0] == '\0' || strcmp(chrom, wi->chrom) < 0) 
			finish = -1;

		if (wi->finish <= finish)
			continue;
		else if (wi->start < finish) 
			wi->start = finish;

		return;
	}

	wi->done = true;
}

static void BedFileReaderOpenFile(BigWiggleReaderData * data, char * f) {
	data->bwf = bigBedFileOpen(f);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;
	data->uncompress = (data->bwf->uncompressBufSize > 0);

	data->chrom = data->chromList = bbiChromList(data->bwf);
	BigBedReaderEnterChromosome(data);
}

void BigBedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigWiggleReaderData * data = (BigWiggleReaderData*) wi->data;
	data->block = data->blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, chrom, start, finish, NULL);
	data->chrom = NULL;
	BigBedReaderEnterRunOfBlocks(data);
}

WiggleIterator * BigBedReader(char * f) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) calloc(1, sizeof(BigWiggleReaderData));
	BedFileReaderOpenFile(data, f);
	WiggleIterator * new = newWiggleIterator(data, &BigBedReaderPop, &BigBedReaderSeek);
	new->value = 1;
	return new;
}	
