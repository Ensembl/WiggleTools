// Copyright [1999-2017] EMBL-European Bioinformatics Institute
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
#include "multiplexer.h"

//////////////////////////////////////////////////////
// Tee operator
//////////////////////////////////////////////////////

#define BLOCK_LENGTH 10000 
#define MAX_OUT_BLOCKS 2

typedef struct BlockData_st {
	char * chroms[BLOCK_LENGTH];
	int starts[BLOCK_LENGTH];
	int finishes[BLOCK_LENGTH];
	double * values;
	int count;
	int width;
	bool bedGraph;
	struct BlockData_st * next;
} BlockData;

typedef struct TeeMultiplexerData_st {
	FILE * infile;
	FILE * outfile;
	Multiplexer * in;
	BlockData * dataBlocks;
	BlockData * lastBlock;
	int count;
	pthread_t threadID;
	pthread_mutex_t continue_mutex;
	pthread_cond_t continue_cond;
	bool done;
	bool bedGraph;
} TeeMultiplexerData;

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
			for (j = 0; j < *finishPtr - *startPtr; j++) {
				int k;
				double * ptr = valuePtr;
				for (k = 0; k < block->width; k++)
					fprintf(outfile, "\t%lf", *(ptr++));
				fprintf(outfile, "\n");
			}
			valuePtr += block->width;
		} else if (!infile) {
			// Careful bedgraph lines are 0 based
			fprintf(outfile, "%s\t%i\t%i", *chromPtr, *startPtr-1, *finishPtr-1);
			int k;
			for (k = 0; k < block->width; k++)
				fprintf(outfile, "\t%lf", *(valuePtr++));
			fprintf(outfile, "\n");
		} else {
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
			int k;
			for (k = strlen(buffer)-1; k >= 0; k--) {
				if (buffer[k] == '\n' || buffer[k] == '\r')
					buffer[k] = '\0';
				else
					break;
			}

			// Print out
			fprintf(outfile, "%s", buffer);
			for (k = 0; k < block->width; k++)
				fprintf(outfile, "\t%lf", *(valuePtr++));
			fprintf(outfile, "\n");
		}

		lastChrom = *chromPtr;
		lastFinish = *finishPtr;
		chromPtr++;
		startPtr++;
		finishPtr++;
	}
}

static bool goToNextBlock(TeeMultiplexerData * data) {
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
	free(ptr->values);
	free(ptr);
	return false;
}

static void * printToFile(void * args) {
	TeeMultiplexerData * data = (TeeMultiplexerData *) args;

	// Wait for first block to arrive
	pthread_mutex_lock(&data->continue_mutex);
	if (data->count == 0 && !data->done) 
		pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
	pthread_mutex_unlock(&data->continue_mutex);

	if (data->count < 0)
		return NULL;

	while(data->dataBlocks) {
		printBlock(data->infile, data->outfile, data->dataBlocks);
		if (goToNextBlock(data))
			return NULL;
	}
	return NULL;
}

static void TeeMultiplexerPop(Multiplexer * multi) {
	TeeMultiplexerData * data = (TeeMultiplexerData *) multi->data;
	Multiplexer * in = data->in;
	if (!data->in->done) {
		multi->chrom = in->chrom;
		multi->start = in->start;
		multi->finish = in->finish;
		// No need to copy values, pointer points to the source multipliexers values' array
		//multi->value = in->value;

		if (data->threadID) {
			int index = data->lastBlock->count;
			data->lastBlock->chroms[index] =  in->chrom;
			data->lastBlock->starts[index] =  in->start;
			data->lastBlock->finishes[index] =  in->finish;
			int i;
			double * ptr = data->lastBlock->values + (index * multi->count);
			for (i = 0; i < multi->count; i++)
				*(ptr++) = in->values[i];

			if (++data->lastBlock->count >= BLOCK_LENGTH) {
				// Communications
				pthread_mutex_lock(&data->continue_mutex);
				data->count++;
				pthread_cond_signal(&data->continue_cond);
				if (data->count > MAX_OUT_BLOCKS)
					pthread_cond_wait(&data->continue_cond, &data->continue_mutex);
				pthread_mutex_unlock(&data->continue_mutex);

				data->lastBlock->next = (BlockData*) calloc(1, sizeof(BlockData));
				data->lastBlock->next->values = (double*) calloc(BLOCK_LENGTH * multi->count, sizeof(double));
				data->lastBlock->next->width = in->count;
				data->lastBlock = data->lastBlock->next;
				data->lastBlock->bedGraph = data->bedGraph;
			}
		}
		popMultiplexer(in);
	} else if (data->threadID) {
		pthread_mutex_lock(&data->continue_mutex);
		data->count++;
		data->done = true;
		pthread_cond_signal(&data->continue_cond);
		pthread_mutex_unlock(&data->continue_mutex);
		multi->done = true;
		pthread_join(data->threadID, NULL);
	}
}

static void launchWriter(TeeMultiplexerData * data, int width) {
	// Initialize variables
	data->count = 0;
	data->done = false;
	pthread_cond_init(&data->continue_cond, NULL);
	pthread_mutex_init(&data->continue_mutex, NULL);
	data->dataBlocks = data->lastBlock = (BlockData*) calloc(1, sizeof(BlockData));
	data->dataBlocks->values = (double*) calloc(BLOCK_LENGTH * width, sizeof(double));
	data->dataBlocks->width = width;
	data->lastBlock->bedGraph = data->bedGraph;

	// Launch pthread
	int err = pthread_create(&data->threadID, NULL, &printToFile, data);
	if (err) {
		fprintf(stderr, "Could not create new thread %i\n", err);
		exit(1);
	}
}

static void killWriter(TeeMultiplexerData * data) {
	BlockData * block;
	
	if (!data->threadID)
		return;

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
		free(block->values);
		free(block);
	}

	data->dataBlocks = NULL;
	data->lastBlock = NULL;
}

static void TeeMultiplexerSeek(Multiplexer * multi, const char * chrom, int start, int finish) {
	TeeMultiplexerData * data = (TeeMultiplexerData *) multi->data;
	killWriter(data);
	fflush(data->outfile);
	seekMultiplexer(data->in, chrom, start, finish);
	multi->done = false;
	launchWriter(data, multi->count);
	popMultiplexer(multi);
}

Multiplexer * TeeMultiplexer(Multiplexer * in, FILE * outfile, bool bedGraph, bool holdFire) {
	TeeMultiplexerData * data = (TeeMultiplexerData *) calloc(1, sizeof(TeeMultiplexerData));
	data->in = in;
	data->outfile = outfile;
	data->bedGraph = bedGraph;
	// Hold fire means that you wait for the first seek before doing any writing
	if (!holdFire)
		launchWriter(data, in->count);

	Multiplexer * res = newCoreMultiplexer(data, in->count, &TeeMultiplexerPop, &TeeMultiplexerSeek);
	res->values = in->values;
	res->inplay = in->inplay;
	res->default_values = in->default_values;
	popMultiplexer(res);
	return res;
}

void toStdoutMultiplexer(Multiplexer * in, bool bedGraph, bool holdFire) {
	runMultiplexer(TeeMultiplexer(in, stdout, bedGraph, holdFire));
}

//////////////////////////////////////////////////////////
// Paste Iterator
//////////////////////////////////////////////////////////

Multiplexer * PasteMultiplexer(Multiplexer * in, FILE * infile, FILE * outfile, bool holdFire) {
	TeeMultiplexerData * data = (TeeMultiplexerData *) calloc(1, sizeof(TeeMultiplexerData));
	data->in = in;
	data->infile = infile;
	data->bedGraph = true;
	data->outfile = outfile;
	// Hold fire means that you wait for the first seek before doing any writing
	if (!holdFire)
		launchWriter(data, in->count);

	Multiplexer * res = newCoreMultiplexer(data, in->count, &TeeMultiplexerPop, &TeeMultiplexerSeek);
	res->values = in->values;
	res->default_values = in->default_values;
	res->inplay = in->inplay;
	popMultiplexer(res);
	return res;
}
