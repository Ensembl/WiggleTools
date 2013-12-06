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

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include <stdint.h>

#include "wiggleIterator.h"

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

static bool appendNewBlock(BinaryFileReaderData * data, BlockData * block) {
	pthread_mutex_lock(&data->count_mutex);
	if (data->count < 0) {
		// Kill signal received
		pthread_mutex_unlock(&data->count_mutex);
		return true;
	}
	// Add block to the end of the list
	if (data->lastBlock) 
		data->lastBlock->next = block;
	else
		data->blocks = block;
	data->lastBlock = block;
	// Signal to reader
	data->count++;
	pthread_cond_signal(&data->count_cond);
	// Pause if too many blocks were read
	if (data->count > MAX_BLOCKS) 
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	pthread_mutex_unlock(&data->count_mutex);
	return false;
}

static void mustRead(void * ptr, size_t size, size_t nelem, FILE * file) {
	if (fread(ptr, size, nelem, file) < nelem) {
		fprintf(stderr, "Incomplete binary file\n");
		exit(1);
	}
}

static bool readNextBlock(BinaryFileReaderData * data, char ** lastChrom, int * lastStart, bool * pointByPoint) {
	BlockData * block = calloc(1, sizeof(BlockData));
	char ** chromPtr = block->chroms;
	int * startPtr = block->starts;
	int * finishPtr = block->finishes;
	double * valuePtr = block->values;
	char c;
	bool startSet = true;
	int32_t holder;
	float holder2;

	for (block->count = 0; block->count < BLOCK_LENGTH; block->count++) {
		// Check whether file finished
		if (fread(&c, 1, 1, data->file) == 0) {
			appendNewBlock(data, block);
			// Increment counter to push the reader into NULL block
			pthread_mutex_lock(&data->count_mutex);
			data->count++;
			pthread_mutex_unlock(&data->count_mutex);
			return true;
		}
		
		if (*pointByPoint)
			startSet = false;

		// Read header
		if (c == 1) {
			// Point by Point
			*pointByPoint = true;
			mustRead(&c, 1, 1, data->file);
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
			mustRead(&holder, sizeof(int32_t), 1, data->file);
			*startPtr = holder;
			startSet = true;
		} else if (c == 2) {
			// Normal
			*pointByPoint = false;
			mustRead(&c, 1, 1, data->file);
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
		}

		// Set chromosome if unset in header
		if (*chromPtr == NULL) 
			*chromPtr = *lastChrom;

		// Read coords
		if (*pointByPoint) {
			if (!startSet)
				*startPtr = *lastStart + 1;
			*finishPtr = *startPtr + 1;
		} else {
			mustRead(&holder, sizeof(int32_t), 1, data->file);
			*startPtr = holder;
			mustRead(&holder, sizeof(int32_t), 1, data->file);
			*finishPtr = holder;
		}

		// Read value
		mustRead(&holder2, sizeof(float), 1, data->file);
		*valuePtr = holder2;

		// Record stuff
		*lastChrom = *chromPtr;
		*lastStart = *startPtr;

		// Step ahead
		chromPtr++;
		startPtr++;
		finishPtr++;
		valuePtr++;
	}
	
	return appendNewBlock(data, block);
}

static void * readBinaryFile(void * args) {
	BinaryFileReaderData * data = (BinaryFileReaderData *) args;
	char * lastChrom = NULL;
	int lastStart = -1;
	bool pointByPoint = false;

	while (true)
		if (readNextBlock(data, &lastChrom, &lastStart, &pointByPoint))
			return NULL;

	return NULL;
}

static void launchDownloader(BinaryFileReaderData * data) {
	pthread_mutex_init(&data->count_mutex, NULL);
	pthread_cond_init(&data->count_cond, NULL);

	int err = pthread_create(&data->threadID, NULL, &readBinaryFile, data);
	if (err) {
		fprintf(stderr, "Could not create new thread %i\n", err);
		exit(1);
	}
}

static void destroyBlockData(BlockData * block) {
	//int pos;

	//for (pos = 0; pos < block->count; pos++)
	//	free(block->chroms[pos]);
	free(block);
}

static void killDownloader(BinaryFileReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	data->count = -1;
	pthread_mutex_unlock(&data->count_mutex);
	pthread_join(data->threadID, NULL);

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
	static int i = 0;
	i++;

	pthread_mutex_lock(&data->count_mutex);
	if (--data->count == 0) 
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	data->blocks = data->blocks->next;
	if (!data->blocks)
		data->lastBlock = NULL;
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
		killDownloader(data);
		wi->done = true;
		return;
	}

        if (data->index == data->blocks->count)
	    BinaryFileReaderGoToNextBlock(data);

	if (!data->blocks) {
		killDownloader(data);
		wi->done = true;
		return;
	}
		
	wi->chrom = data->blocks->chroms[data->index];
	wi->start = data->blocks->starts[data->index];
	wi->finish = data->blocks->finishes[data->index];
	wi->value = data->blocks->values[data->index];

	data->index++;

	if (data->stop > 0 && (wi->start > data->stop)) {
		killDownloader(data);
		wi->done = true;
		return;
	}

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
	if (!data->file) {
		fprintf(stderr, "Could not open %s\n", f);
		exit(1);
	}
	launchDownloader(data);
	BinaryFileReaderEnterBlock(data);

	return newWiggleIterator(data, &BinaryFileReaderPop, &BinaryFileReaderSeek);
}	
