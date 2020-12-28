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

#include <string.h>
#include "bufferedReader.h"

static int MAX_HEAD_START = 3;
static int BLOCK_SIZE = 10000;

typedef struct blockData_st {
	char **chrom;
	int * start;
	int * finish;
	double * value;
	int count;
	struct blockData_st * next;
} BlockData;

struct bufferedReaderData_st {
	pthread_t downloaderThreadID;
	BlockData * blockData, *lastBlockData;
	pthread_mutex_t count_mutex;
	pthread_cond_t count_cond;
	int blockCount;
	int readIndex;
	void * readerData;
	bool killed;
};

static bool declareNewBlock(BufferedReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);

	if (data->blockCount > MAX_HEAD_START) {
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	} 
	if (data->blockCount < 0) {
		pthread_mutex_unlock(&data->count_mutex);
		return true;
	}
	data->blockCount++;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
	return false;
}

static BlockData * createBlockData() {
	BlockData * new = (BlockData * ) calloc(1, sizeof(BlockData));
	new->chrom = (char **) calloc(BLOCK_SIZE, sizeof(char*));
	new->start = (int *) calloc(BLOCK_SIZE, sizeof(int));
	new->finish = (int *) calloc(BLOCK_SIZE, sizeof(int));
	new->value = (double *) calloc(BLOCK_SIZE, sizeof(double));
	return new;
}

bool pushValuesToBuffer(BufferedReaderData * data, const char * chrom, int start, int finish, double value) {

	if (data->blockData == NULL)
		data->lastBlockData = data->blockData = createBlockData();
	else if (data->lastBlockData->count == BLOCK_SIZE) {
		data->lastBlockData->next = createBlockData();
		data->lastBlockData = data->lastBlockData->next;
		if (declareNewBlock(data))
			return true;
	}

	int index = data->lastBlockData->count;
	data->lastBlockData->chrom[index] = chrom;
	data->lastBlockData->start[index] = start;
	data->lastBlockData->finish[index] = finish;
	data->lastBlockData->value[index] = value;
	data->lastBlockData->count++;
	return false;
}

void endBufferedSignal(BufferedReaderData * data) {
	declareNewBlock(data);
	declareNewBlock(data);
}

static void destroyBlockData(BlockData * data) {
	free(data->chrom);
	free(data->start);
	free(data->finish);
	free(data->value);
	free(data);
}

static void waitForNextBlock(BufferedReaderData * data) {
	pthread_mutex_lock(&data->count_mutex);
	// Check whether allowed to step forward
	if (data->blockCount == 0) {
		pthread_cond_wait(&data->count_cond, &data->count_mutex);
	}
	// Signal freed memory
	data->blockCount--;
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
}

static void goToNextBlock(BufferedReaderData * data) {
	BlockData * prevBlockData = data->blockData;
	data->blockData = data->blockData->next;
	data->readIndex = 0;
	destroyBlockData(prevBlockData);
}

void launchBufferedReader(void * (* readFileFunction)(void *), void * f_data, BufferedReaderData ** buf_data) {
	BufferedReaderData * data = calloc(1, sizeof(BufferedReaderData));
	*buf_data = data;
	data->readIndex = 0;
	data->readerData = f_data;

	pthread_mutex_init(&data->count_mutex, NULL);
	pthread_cond_init(&data->count_cond, NULL);

	int err = pthread_create(&(data->downloaderThreadID), NULL, readFileFunction, f_data);
	if (err) {
		fprintf(stderr, "Could not create new thread %i\n", err);
		abort();
	}

	waitForNextBlock(data);
}

void killBufferedReader(BufferedReaderData * data) {
	if (data->killed)
		return;

	pthread_mutex_lock(&data->count_mutex);
	data->blockCount = -1;
	// Send a signal in case the slave is waiting somewhere
	pthread_cond_signal(&data->count_cond);
	pthread_mutex_unlock(&data->count_mutex);
	pthread_join(data->downloaderThreadID, NULL);

	pthread_mutex_destroy(&data->count_mutex);
	pthread_cond_destroy(&data->count_cond);

	while (data->blockData) {
		BlockData * prevData = data->blockData;
		data->blockData = data->blockData->next;
		destroyBlockData(prevData);
	}

	data->lastBlockData = NULL;
	data->blockCount = 0;
	data->killed = true;
}

void BufferedReaderPop(WiggleIterator * wi, BufferedReaderData * data) {
	if (wi->done)
		return;
	else if (data == NULL || data->blockData == NULL) {
		wi->done = true;
		return;
	} else if (data->readIndex == data->blockData->count) {
		waitForNextBlock(data);
		goToNextBlock(data);
		if (data->blockData == NULL) {
			killBufferedReader(data);
			wi->done = true;
			return;
		}
	} 

	int index = data->readIndex;
	wi->chrom = data->blockData->chrom[index];
	wi->start = data->blockData->start[index];
	wi->finish = data->blockData->finish[index];
	wi->value = (double) data->blockData->value[index];
	data->readIndex++;
}


int compare_chrom_lengths(const void * A, const void * B) {
	Chrom_length * cl_A = (Chrom_length *) A;
	Chrom_length * cl_B = (Chrom_length *) B;
	return strcmp(cl_A->chrom, cl_B->chrom);
}
