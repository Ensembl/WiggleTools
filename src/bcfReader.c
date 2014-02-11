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

#include <string.h>
#include <zlib.h>
#include <tabix.h>
#include "wiggleIterator.h"
#include "bufferedReader.h"

#define BUFF_LENGTH 1000


typedef struct bamFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int start, stop;
	BufferedReaderData * bufferedReaderData;
	char buffer[BUFF_LENGTH];

	// Tabix stuff
	tabix_t * tabix_file;
	ti_iter_t tabix_iterator;

	// Gzip file
	gzFile gz_file;
} VCFReaderData;

static char * nextLine(VCFReaderData * data) {
	if (data->tabix_iterator)
		return ti_read(data->tabix_file, data->tabix_iterator, 0);
	else
		return gzgets(data->gz_file, data->buffer, BUFF_LENGTH);
}

static void * downloadTabixFile(void * args) {
	VCFReaderData * data = (VCFReaderData *) args;
	char * line;
	char * last_chrom = "";

	while ((line = nextLine(data))) {
		if (line[0] == '#')
			continue;

		char * chrom = strtok(line, "\t");
		if (strcmp(chrom, last_chrom)) {
			last_chrom = calloc(strlen(chrom) + 1, sizeof(char));
			strcpy(last_chrom, chrom);
		}
		int pos = atoi(strtok(NULL, "\t"));

		pushValuesToBuffer(data->bufferedReaderData, last_chrom, pos, pos+1, 1);
		if (data->tabix_iterator)
			free(line);
	}

	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void OpenTabixFile(VCFReaderData * data, char * filename) {
	data->tabix_file = ti_open(filename, NULL);
	data->gz_file = gzopen(filename, "r");
}

void killVCFReader(VCFReaderData * data) {
	if (data->tabix_iterator) 
		ti_iter_destroy(data->tabix_iterator);
}

void closeTabixFile(VCFReaderData * data) {
	ti_close(data->tabix_file);
	gzclose(data->gz_file);
}

void VCFReaderPop(WiggleIterator * wi) {
	VCFReaderData * data = (VCFReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void BcfReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	killBufferedReader(wi->data);

	VCFReaderData * data = (VCFReaderData *) wi->data;
	data->tabix_iterator = ti_query(data->tabix_file, chrom, start, finish);
	launchBufferedReader(&downloadTabixFile, &killVCFReader, data, &(data->bufferedReaderData));
	wi->done = false;
	VCFReaderPop(wi);

	while (strcmp(wi->chrom, chrom) < 0 || wi->finish <= start)
		VCFReaderPop(wi);

	data->chrom = chrom;
	data->stop = finish;
}

WiggleIterator * BcfReader(char * filename) {
	VCFReaderData * data = (VCFReaderData *) calloc(1, sizeof(VCFReaderData));
	OpenTabixFile(data, filename);
	launchBufferedReader(&downloadTabixFile, &killVCFReader, data, &(data->bufferedReaderData));
	return newWiggleIterator(data, &VCFReaderPop, &BcfReaderSeek);
}
