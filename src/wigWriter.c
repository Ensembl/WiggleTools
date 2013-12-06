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
#include <string.h>
#include <pthread.h>
#include <stdint.h>

// Local header
#include "wiggleIterator.h"

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
	bool bedGraph;
	struct BlockData_st * next;
} BlockData;

typedef struct TeeWiggleIteratorData_st {
	FILE * infile;
	FILE * outfile;
	WiggleIterator * iter;
	BlockData * dataBlocks;
	BlockData * lastBlock;
	int count;
	pthread_t threadID;
	pthread_mutex_t continue_mutex;
	pthread_cond_t continue_cond;
	bool done;
	bool binary;
	bool bedGraph;
} TeeWiggleIteratorData;

static void printBlock(FILE * infile, FILE * outfile, BlockData * block) {
	int i, j;
	bool pointByPoint = false;
	bool makeHeader=false;
	char ** chromPtr = block->chroms;
	int * startPtr = block->starts;
	int * finishPtr = block->finishes;
	double * valuePtr = block->values;
	char * lastChrom = NULL;
	int lastFinish = -1;
	char buffer[5000];

	for (i = 0; i < block->count; i++) {
		// Change mode
		if (!block->bedGraph && *finishPtr - *startPtr < 2 && !pointByPoint) {
			pointByPoint = true;
			makeHeader = true;
		} else if (*finishPtr - *startPtr > 5 && pointByPoint) {
			pointByPoint = false;
		}

		if (pointByPoint) {
			if (makeHeader || (pointByPoint && (lastChrom != *chromPtr || *startPtr > lastFinish)))
				fprintf(outfile, "fixedStep chrom=%s start=%i step=1\n", *chromPtr, *startPtr);
			makeHeader = false;
			for (j = 0; j < *finishPtr - *startPtr; j++)
				fprintf(outfile, "%lf\n", *valuePtr);
		} else if (!infile)
			// Careful bedgraph lines are 0 based
			fprintf(outfile, "%s\t%i\t%i\t%lf\n", *chromPtr, *startPtr-1, *finishPtr-1, *valuePtr);
		else {
			// Read next line in infile
			if (!fgets(buffer, 5000, infile)) {
				fprintf(stderr, "Could not paste data to file lines, inconsistent number of lines.\n");
				exit(1);
			}

			// Skip empty lines and metadata lines:
			while (! (strlen(buffer) && strncmp(buffer, "track", 5) && strncmp(buffer, "browser", 7))) {
				if (!fgets(buffer, 5000, infile)) {
					fprintf(stderr, "Could not paste data to file lines, inconsistent number of lines.\n");
					exit(1);
				}
			}

			// Strip end of line symbols
			int i;
			for (i = strlen(buffer)-1; i >= 0; i--) {
				if (buffer[i] == '\n' || buffer[i] == '\r')
					buffer[i] = '\0';
				else
					break;
			}
			// Print out
			fprintf(outfile, "%s\t%lf\n", buffer, *valuePtr);
		}

		lastChrom = *chromPtr;
		lastFinish = *finishPtr;
		chromPtr++;
		startPtr++;
		finishPtr++;
		valuePtr++;
	}
}

static void printBinaryBlock(FILE * infile, FILE * outfile, BlockData * block) {
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
				fwrite(&flag, sizeof(char), 1, outfile);
				if (*chromPtr == lastChrom) {
					fwrite(emptyString, sizeof(char), 1, outfile);
				} else {
					fwrite(*chromPtr, sizeof(char), strlen(*chromPtr) + 1, outfile);
				}
				holder = (int32_t) *startPtr;
				fwrite(&holder, sizeof(int32_t), 1, outfile);
			} else {
				flag = (char) 2;
				fwrite(&flag, sizeof(char), 1, outfile);
				if (*chromPtr == lastChrom) {
					fwrite(emptyString, sizeof(char), 1, outfile);
				} else
					fwrite(*chromPtr, sizeof(char), strlen(*chromPtr) + 1, outfile);
			}
		} else {
			fwrite(emptyString, sizeof(char), 1, outfile);
		}

		// Data
		if (!pointByPoint) {
			holder = (int32_t) *startPtr;
			fwrite(&holder, sizeof(int32_t), 1, outfile);
			holder = (int32_t) *finishPtr;
			fwrite(&holder, sizeof(int32_t), 1, outfile);
			holder2 = *valuePtr;
			fwrite(&holder2, sizeof(float), 1, outfile);
		} else {
			int j;
			holder2 = *valuePtr;
			fwrite(&holder2, sizeof(float), 1, outfile);

			for (j = 1; j < *finishPtr - *startPtr; j++) {
				fwrite(emptyString, sizeof(char), 1, outfile);
				holder2 = *valuePtr;
				fwrite(&holder2, sizeof(float), 1, outfile);
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
	BlockData * ptr = data->dataBlocks;
	static int i = 0;
	i++;

	pthread_mutex_lock(&data->continue_mutex);
	// Received kill signal
	if (data->count < 0)
		return true;

	// Check that there is work left
	data->count--;
	if (data->count == 0 && !data->done) 
		pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
	pthread_cond_signal(&data->continue_cond);
	pthread_mutex_unlock(&data->continue_mutex);

	// Step forward
	data->dataBlocks = data->dataBlocks->next;
	free(ptr);
	return false;
}

static void * printToFile(void * args) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) args;

	// Wait for first block to arrive
	pthread_mutex_lock(&data->continue_mutex);
	if (data->count == 0 && !data->done) 
		pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
	pthread_mutex_unlock(&data->continue_mutex);

	if (data->count < 0)
		return NULL;

	while(data->dataBlocks) {
		if (data->binary)
			printBinaryBlock(data->infile, data->outfile, data->dataBlocks);
		else
			printBlock(data->infile, data->outfile, data->dataBlocks);
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

		int index = data->lastBlock->count;
		data->lastBlock->chroms[index] =  iter->chrom;
		data->lastBlock->starts[index] =  iter->start;
		data->lastBlock->finishes[index] =  iter->finish;
		data->lastBlock->values[index] =  iter->value;
		if (++data->lastBlock->count >= BLOCK_LENGTH) {
			// Communications
			pthread_mutex_lock(&data->continue_mutex);
			data->count++;
			pthread_cond_signal(&data->continue_cond);
			if (data->count > MAX_OUT_BLOCKS)
				pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
			pthread_mutex_unlock(&data->continue_mutex);

			data->lastBlock->next = (BlockData*) calloc(1, sizeof(BlockData));
			data->lastBlock = data->lastBlock->next;
			data->lastBlock->bedGraph = data->bedGraph;
		}
		pop(iter);
	} else {
		pthread_mutex_lock(&data->continue_mutex);
		data->count++;
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
	data->dataBlocks = data->lastBlock = (BlockData*) calloc(1, sizeof(BlockData));
	data->lastBlock->bedGraph = data->bedGraph;

	// Launch pthread
	int err = pthread_create(&data->threadID, NULL, &printToFile, data);
	if (err) {
		fprintf(stderr, "Could not create new thread %i\n", err);
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

	while (data->dataBlocks) {
		block = data->dataBlocks;
		data->dataBlocks = block->next;
		free(block);
	}

	data->dataBlocks = NULL;
	data->lastBlock = NULL;
}

void TeeWiggleIteratorSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) wi->data;
	killWriter(data);
	fflush(data->outfile);
	fseek(data->outfile, 0, SEEK_SET);
	seek(data->iter, chrom, start, finish);
	wi->done = false;
	launchWriter(data);
	pop(wi);
}

WiggleIterator * BinaryTeeWiggleIterator(WiggleIterator * i, FILE * outfile, bool bedGraph) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) calloc(1, sizeof(TeeWiggleIteratorData));
	data->iter = CompressionWiggleIterator(i);
	data->outfile = outfile;
	data->binary = true;
	data->bedGraph = bedGraph;
	launchWriter(data);

	return newWiggleIterator(data, &TeeWiggleIteratorPop, &TeeWiggleIteratorSeek);
}

WiggleIterator * TeeWiggleIterator(WiggleIterator * i, FILE * outfile, bool bedGraph) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) calloc(1, sizeof(TeeWiggleIteratorData));
	data->iter = CompressionWiggleIterator(i);
	data->outfile = outfile;
	data->bedGraph = bedGraph;
	launchWriter(data);

	return newWiggleIterator(data, &TeeWiggleIteratorPop, &TeeWiggleIteratorSeek);
}

void toBinaryFile(WiggleIterator * wi, char * filename, bool bedGraph) {
	FILE * file = fopen(filename, "wb");
	if (!file) {
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(1);
	}
	runWiggleIterator(BinaryTeeWiggleIterator(wi, file, bedGraph));
}

void toFile(WiggleIterator * wi, char * filename, bool bedGraph) {
	FILE * file = fopen(filename, "w");
	if (!file) {
		fprintf(stderr, "Could not open file %s\n", filename);
		exit(1);
	}
	runWiggleIterator(TeeWiggleIterator(wi, file, bedGraph));
}

void toStdout(WiggleIterator * wi, bool bedGraph) {
	runWiggleIterator(TeeWiggleIterator(wi, stdout, bedGraph));
}

//////////////////////////////////////////////////////////
// Paste Iterator
//////////////////////////////////////////////////////////

WiggleIterator * PasteWiggleIterator(WiggleIterator * i, FILE * infile, FILE * outfile) {
	TeeWiggleIteratorData * data = (TeeWiggleIteratorData *) calloc(1, sizeof(TeeWiggleIteratorData));
	data->iter = i;
	data->infile = infile;
	data->bedGraph = true;
	data->outfile = outfile;
	launchWriter(data);

	return newWiggleIterator(data, &TeeWiggleIteratorPop, NULL);
}
