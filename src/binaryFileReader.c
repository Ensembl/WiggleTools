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
#include <pthread.h>
#include <string.h>

#include "wiggleIterators.h"

static int MAX_BLOCKS = 100;
#define BLOCK_LENGTH 10000

void setBinaryMaxBlocks(int value) {
	MAX_BLOCKS = value;
}

typedef struct BlockData_st {
	char * chroms[BLOCK_LENGTH];
	int starts[BLOCK_LENGTH];
	int finishes[BLOCK_LENGTH];
	double values[BLOCK_LENGTH];
	int count;
	struct BlockData_st * next;
} BlockData;

typedef struct BinaryFileReaderData_st {
	FILE * file;
	BlockData * blocks, *lastBlock;
	int count;
	pthread_t threadID;
	pthread_mutex_t count_mutex;
	pthread_cond_t count_cond;
	char * chrom;
	int stop;
	int index;
} BinaryFileReaderData;

static void appendNewBlock(BinaryFileReaderData * data, BlockData * block) {
	pthread_mutex_lock(&data->count_mutex);
	if (data->lastBlock) 
		data->lastBlock->next = block;
	else
		data->blocks = block;
	data->lastBlock = block;
	pthread_cond_signal(&data->count_cond);
	data->count++;
	if (data->count > MAX_BLOCKS)
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	pthread_mutex_unlock(&data->count_mutex);
}

static void mustRead(void * ptr, size_t size, size_t nelem, FILE * file) {
	if (fread(ptr, size, nelem, file) < nelem) {
		printf("Incomplete binary file\n");
		exit(1);
	}
}

static bool readNextBlock(BinaryFileReaderData * data) {
	BlockData * block = calloc(1, sizeof(BlockData));
	int i;
	char ** chromPtr = block->chroms;
	int * startPtr = block->starts;
	int * finishPtr = block->finishes;
	double * valuePtr = block->values;
	char c;

	for (i = 0; i < BLOCK_LENGTH; i++) {
		if (fread(&c, 1, 1, data->file) == 0) {
			block->count = i;
			appendNewBlock(data, block);
			return true;
		}

		if (c) {
			*chromPtr = calloc(1000, 1);
			char * ptr = * chromPtr;
			int pos;
			for (pos = 0; pos < 1000 && c; pos++) {
				*ptr = c;
				mustRead(&c, 1, 1, data->file);
				ptr++;
			}
			
			while (c)
				mustRead(&c, 1, 1, data->file);
		} 

		mustRead(startPtr, sizeof(int), 1, data->file);
		mustRead(finishPtr, sizeof(int), 1, data->file);
		mustRead(valuePtr, sizeof(double), 1, data->file);

		chromPtr++;
		startPtr++;
		finishPtr++;
		valuePtr++;
	}
	block->count = i;
	
	appendNewBlock(data, block);
	return false;

}

static void * readBinaryFile(void * args) {
	BinaryFileReaderData * data = (BinaryFileReaderData *) args;

	while (true)
		if (readNextBlock(data))
			return NULL;

	return NULL;
}

static void launchDownloader(BinaryFileReaderData * data) {
	pthread_mutex_init(&data->count_mutex, NULL);
	pthread_cond_init(&data->count_cond, NULL);

	int err = pthread_create(&data->threadID, NULL, &readBinaryFile, data);
	if (err) {
		printf("Could not create new thread %i\n", err);
		abort();
	}
	pthread_detach(data->threadID);
}

static void destroyBlockData(BlockData * block) {
	//int pos;

	//for (pos = 0; pos < block->count; pos++)
	//	free(block->chroms[pos]);
	free(block);
}

static void killDownloader(BinaryFileReaderData * data) {
	if (data->threadID)
		pthread_cancel(data->threadID);

	pthread_mutex_destroy(&data->count_mutex);
	pthread_cond_destroy(&data->count_cond);

	while (data->blocks) {
		BlockData * prevData = data->blocks;
		data->blocks = data->blocks->next;
		destroyBlockData(prevData);
	}

	data->count = 0;
}

static void BinaryFileReaderEnterBlock(BinaryFileReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	if (data->count == 0)
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	pthread_mutex_unlock(&data->count_mutex);

	data->index = 0;
}

static void BinaryFileReaderGoToNextBlock(BinaryFileReaderData * data) {
	BlockData * ptr = data->blocks;

	pthread_mutex_lock(&data->count_mutex);
	if (data->count == 0)
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	data->blocks = data->blocks->next;
	if (!data->blocks)
		data->lastBlock = NULL;
	data->count--;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);

	data->index = 0;
	destroyBlockData(ptr);
}

void BinaryFileReaderPop(WiggleIterator * wi) {
	BinaryFileReaderData * data;

	if (wi->done)
		return;

	data = (BinaryFileReaderData*) wi->data;

	if (!data->blocks) {
		wi->done = true;
		return;
	}
		
	if (data->blocks->chroms[data->index])
		wi->chrom = data->blocks->chroms[data->index];
	wi->start = data->blocks->starts[data->index];
	wi->finish = data->blocks->finishes[data->index];
	wi->value = data->blocks->values[data->index];

	if (data->stop > 0 && (wi->start > data->stop)) {
		wi->done = true;
		return;
	}

        if (++(data->index) == data->blocks->count)
	    BinaryFileReaderGoToNextBlock(data);

}

void BinaryFileReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BinaryFileReaderData * data = (BinaryFileReaderData *) wi->data; 

	killDownloader(data);
	data->chrom = chrom;
	data->stop = finish;
	launchDownloader(data);
	BinaryFileReaderEnterBlock(data);
	wi->done = false;
	BinaryFileReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BinaryFileReaderPop(wi);
}

WiggleIterator * BinaryFileReader(char * f) {
	BinaryFileReaderData * data = (BinaryFileReaderData *) calloc(1, sizeof(BinaryFileReaderData));
	data->file = fopen(f, "rb");
	launchDownloader(data);
	BinaryFileReaderEnterBlock(data);

	return newWiggleIterator(data, &BinaryFileReaderPop, &BinaryFileReaderSeek);
}	
