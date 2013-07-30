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
#include <pthread.h>
#include <stdint.h>

// Local header
#include "wiggleTools.h"
#include "wiggleIterators.h"

//////////////////////////////////////////////////////
// Tee operator
//////////////////////////////////////////////////////

#define BLOCK_LENGTH 10000 
#define MAX_OUT_BLOCKS 2

typedef struct BlockData_st {
	char * chroms[BLOCK_LENGTH];
	int starts[BLOCK_LENGTH];
	int finishes[BLOCK_LENGTH];
	double values[BLOCK_LENGTH];
	int count;
	struct BlockData_st * next;
} BlockData;

typedef struct TeeWiggleIteratorData_st {
	FILE * file;
	WiggleIterator * iter;
	BlockData * finishedBlocks;
	BlockData * lastBlock;
	BlockData * fillingBlock;
	int count;
	pthread_t threadID;
	pthread_mutex_t continue_mutex;
	pthread_cond_t continue_cond;
	bool done;
	bool binary;
} TeeWiggleIteratorData;

static void appendFinishedBlock(TeeWiggleIteratorData * data) {
	if (data->lastBlock)
		data->lastBlock->next = data->fillingBlock;
	else
		data->finishedBlocks = data->fillingBlock;
	data->lastBlock = data->fillingBlock;
	data->count++;
}

static void printBlock(FILE * file, BlockData * block) {
	int i, j;
	bool pointByPoint = false;
	bool makeHeader=false;
	char ** chromPtr = block->chroms;
	int * startPtr = block->starts;
	int * finishPtr = block->finishes;
	double * valuePtr = block->values;
	char * lastChrom = NULL;
	int lastFinish = -1;

	for (i = 0; i < block->count; i++) {
		// Change mode
		if (*finishPtr - *startPtr < 2 && !pointByPoint) {
			pointByPoint = true;
			makeHeader = true;
		} else if (*finishPtr - *startPtr > 5 && pointByPoint) {
			pointByPoint = false;
		}

		if (pointByPoint) {
			if (makeHeader || (pointByPoint && (lastChrom != *chromPtr || *startPtr > lastFinish)))
				fprintf(file, "fixedStep chrom=%s start=%i step=1\n", *chromPtr, *startPtr);
			makeHeader = false;
			for (j = 0; j < *finishPtr - *startPtr; j++)
				fprintf(file, "%lf\n", *valuePtr);
		} else
			fprintf(file, "%s\t%i\t%i\t%lf\n", *chromPtr, *startPtr-1, *finishPtr-1, *valuePtr);

		lastChrom = *chromPtr;
		lastFinish = *finishPtr;
		chromPtr++;
		startPtr++;
		finishPtr++;
		valuePtr++;
	}
}

static void printBinaryBlock(FILE * file, BlockData * block) {
	int i;
	char ** chromPtr = block->chroms;
	int * startPtr = block->starts;
	int * finishPtr = block->finishes;
	double * valuePtr = block->values;
	char * emptyString = "";
	char * lastChrom = NULL;
	int lastFinish = -1;
	bool pointByPoint = false;
	bool makeHeader;
	char flag;
	int32_t holder;
	float holder2;

	for (i = 0; i < block->count; i++) {
		// Change mode
		if (*finishPtr - *startPtr < 2 && !pointByPoint) {
			pointByPoint = true;
			makeHeader = true;
		} else if (*finishPtr - *startPtr > 5 && pointByPoint) {
			pointByPoint = false;
			makeHeader = true;
		} else 
			makeHeader = false;

		// Detect discontinuity
		if (lastChrom != *chromPtr || (pointByPoint && *startPtr > lastFinish))
			makeHeader = true;

		// Make header
		if (makeHeader) {
			if (pointByPoint) {
				flag = (char) 1;
				fwrite(&flag, sizeof(char), 1, file);
				if (*chromPtr == lastChrom) {
					fwrite(emptyString, sizeof(char), 1, file);
				} else {
					fwrite(*chromPtr, sizeof(char), strlen(*chromPtr) + 1, file);
				}
				holder = (int32_t) *startPtr;
				fwrite(&holder, sizeof(int32_t), 1, file);
			} else {
				flag = (char) 2;
				fwrite(&flag, sizeof(char), 1, file);
				if (*chromPtr == lastChrom) {
					fwrite(emptyString, sizeof(char), 1, file);
				} else
					fwrite(*chromPtr, sizeof(char), strlen(*chromPtr) + 1, file);
			}
		} else {
			fwrite(emptyString, sizeof(char), 1, file);
		}

		// Data
		if (!pointByPoint) {
			holder = (int32_t) *startPtr;
			fwrite(&holder, sizeof(int32_t), 1, file);
			holder = (int32_t) *finishPtr;
			fwrite(&holder, sizeof(int32_t), 1, file);
			holder2 = *valuePtr;
			fwrite(&holder2, sizeof(float), 1, file);
		} else {
			int j;
			holder2 = *valuePtr;
			fwrite(&holder2, sizeof(float), 1, file);

			for (j = 1; j < *finishPtr - *startPtr; j++) {
				fwrite(emptyString, sizeof(char), 1, file);
				holder2 = *valuePtr;
				fwrite(&holder2, sizeof(float), 1, file);
			}
		}

		lastChrom = *chromPtr;
		lastFinish = *finishPtr;
		chromPtr++;
		startPtr++;
		finishPtr++;
		valuePtr++;
	}
}

static bool goToNextBlock(TeeWiggleIteratorData * data) {
	BlockData * ptr = data->finishedBlocks;
	static int i = 0;
	i++;

	pthread_mutex_lock(&data->continue_mutex);
	if (data->count < 0)
		return true;
	if (!data->finishedBlocks->next && !data->done) 
		pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
	data->finishedBlocks = data->finishedBlocks->next;
	if (!data->finishedBlocks)
		data->lastBlock = NULL;
	data->count--;
	pthread_cond_signal(&data->continue_cond);
	pthread_mutex_unlock(&data->continue_mutex);

	free(ptr);
	return false;
}

static void * printToFile(void * args) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) args;

	// Wait for first block to arrive
	pthread_mutex_lock(&data->continue_mutex);
	if (!data->finishedBlocks && !data->done)
		pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
	pthread_mutex_unlock(&data->continue_mutex);

	while(data->finishedBlocks) {
		if (data->binary)
			printBinaryBlock(data->file, data->finishedBlocks);
		else
			printBlock(data->file, data->finishedBlocks);
		if (goToNextBlock(data))
			return NULL;
	}
	return NULL;
}

void TeeWiggleIteratorPop(WiggleIterator * wi) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) wi->data;
	WiggleIterator * iter = data->iter;
	if (!data->iter->done) {
		wi->chrom = iter->chrom;
		wi->start = iter->start;
		wi->finish = iter->finish;
		wi->value = iter->value;

		int index = data->fillingBlock->count;
		data->fillingBlock->chroms[index] =  iter->chrom;
		data->fillingBlock->starts[index] =  iter->start;
		data->fillingBlock->finishes[index] =  iter->finish;
		data->fillingBlock->values[index] =  iter->value;
		if (++data->fillingBlock->count >= BLOCK_LENGTH) {
			pthread_mutex_lock(&data->continue_mutex);
			appendFinishedBlock(data);
			pthread_cond_signal(&data->continue_cond);
			if (data->count > MAX_OUT_BLOCKS)
				pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
			pthread_mutex_unlock(&data->continue_mutex);
			data->fillingBlock = (BlockData*) calloc(1, sizeof(BlockData));
		}
		pop(iter);
	} else {
		pthread_mutex_lock(&data->continue_mutex);
		appendFinishedBlock(data);
		data->done = true;
		pthread_cond_signal(&data->continue_cond);
		pthread_mutex_unlock(&data->continue_mutex);
		wi->done = true;
		pthread_join(data->threadID, NULL);
	}
}

static void launchWriter(TeeWiggleIteratorData * data) {
	// Initialize variables
	data->count = 0;
	data->done = false;
	pthread_cond_init(&data->continue_cond, NULL);
	pthread_mutex_init(&data->continue_mutex, NULL);
	data->fillingBlock = (BlockData*) calloc(1, sizeof(BlockData));

	// Launch pthread
	int err = pthread_create(&data->threadID, NULL, &printToFile, data);
	if (err) {
		printf("Could not create new thread %i\n", err);
		exit(1);
	}
}

static void killWriter(TeeWiggleIteratorData * data) {
	BlockData * block;

	// Set trap
	pthread_mutex_lock(&data->continue_mutex);
	data->count = -1;
	pthread_cond_signal(&data->continue_cond);
	pthread_mutex_unlock(&data->continue_mutex);

	// Wait for the catch
	pthread_join(data->threadID, NULL);

	// Clear variables
	pthread_cond_destroy(&data->continue_cond);
	pthread_mutex_destroy(&data->continue_mutex);

	while (data->finishedBlocks) {
		block = data->finishedBlocks;
		data->finishedBlocks = block->next;
		free(block);
	}

	data->finishedBlocks = NULL;
	data->lastBlock = NULL;
	
	if (data->fillingBlock) 
		free(data->fillingBlock);

	data->fillingBlock = NULL;
}

void TeeWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) wi->data;
	killWriter(data);
	fflush(data->file);
	fseek(data->file, 0, SEEK_SET);
	seek(data->iter, chrom, start, finish);
	wi->done = false;
	launchWriter(data);
	pop(wi);
}

WiggleIterator * BinaryTeeWiggleIterator(WiggleIterator * i, FILE * file) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) calloc(1, sizeof(TeeWiggleIteratorData));
	data->iter = CompressionWiggleIterator(i);
	data->file = file;
	data->binary = true;
	launchWriter(data);

	return newWiggleIterator(data, &TeeWiggleIteratorPop, &TeeWiggleIteratorSeek);
}

WiggleIterator * TeeWiggleIterator(WiggleIterator * i, FILE * file) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) calloc(1, sizeof(TeeWiggleIteratorData));
	data->iter = CompressionWiggleIterator(i);
	data->file = file;
	launchWriter(data);

	return newWiggleIterator(data, &TeeWiggleIteratorPop, &TeeWiggleIteratorSeek);
}

void toBinaryFile(WiggleIterator * wi, char * filename) {
	FILE * file = fopen(filename, "wb");
	if (!file) {
		printf("Could not open file %s\n", filename);
		exit(1);
	}
	runWiggleIterator(BinaryTeeWiggleIterator(wi, file));
}

void toFile(WiggleIterator * wi, char * filename) {
	FILE * file = fopen(filename, "w");
	if (!file) {
		printf("Could not open file %s\n", filename);
		exit(1);
	}
	runWiggleIterator(TeeWiggleIterator(wi, file));
}

void toStdout(WiggleIterator * wi) {
	runWiggleIterator(TeeWiggleIterator(wi, stdout));
}

