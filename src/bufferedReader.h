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

#ifndef _BUFFERED_READER_H_
#define _BUFFERED_READER_H_

#include <stdlib.h>
#include <pthread.h>
#include "wiggletools.h"
#include "wiggleIterator.h"

typedef struct chrom_length_st {
	char * chrom;
	int length;
} Chrom_length;
typedef struct bufferedReaderData_st BufferedReaderData;

void launchBufferedReader(void * (* readFileFunction)(void *), void * f_data, BufferedReaderData ** buf_data);
bool pushValuesToBuffer(BufferedReaderData * data, const char * chrom, int start, int finish, double value);
void endBufferedSignal(BufferedReaderData * data);
void killBufferedReader(BufferedReaderData * data);
void BufferedReaderPop(WiggleIterator * wi, BufferedReaderData * data);

int compare_chrom_lengths(const void * A, const void * B);
#endif
