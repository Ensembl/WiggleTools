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
#include <zlib.h>
#include "htslib/vcf.h"
#include "wiggleIterator.h"
#include "bufferedReader.h"

#define BUFF_LENGTH 1000


typedef struct bamFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int start, stop;
	BufferedReaderData * bufferedReaderData;

	// BCF stuff
	htsFile * bcf_file;
	bcf_hdr_t * bcf_header;
	hts_itr_t * bcf_iterator;
	hts_idx_t * bcf_index;

	// Gzip file
} BCFReaderData;

static int nextLine(BCFReaderData * data, bcf1_t * holder) {
	if (data->bcf_iterator)
		return bcf_itr_next(data->bcf_file, data->bcf_iterator, holder);
	else
		return bcf_read(data->bcf_file, data->bcf_header, holder);
}

static void * downloadBCFFile(void * args) {
	BCFReaderData * data = (BCFReaderData *) args;
	bcf1_t * vcf_line = bcf_init();

	while (nextLine(data, vcf_line) >= 0) 
		// Note that BCF encoding is 0-based, hence +1s
		if (pushValuesToBuffer(data->bufferedReaderData, bcf_hdr_id2name(data->bcf_header, vcf_line->rid), vcf_line->pos+1, vcf_line->pos+2, 1))
			break;

	endBufferedSignal(data->bufferedReaderData);
	bcf_destroy(vcf_line);
	if (data->bcf_iterator) 
		bcf_itr_destroy(data->bcf_iterator);
	return NULL;
}

void OpenBCFFile(BCFReaderData * data, char * filename) {
	data->bcf_file = bcf_open(filename, "r");
	data->bcf_index = bcf_index_load(filename);
	data->bcf_header = bcf_hdr_read(data->bcf_file);
}

void closeBCFFile(BCFReaderData * data) {
	if (data->bcf_iterator)
		bcf_itr_destroy(data->bcf_iterator);
	bcf_hdr_destroy(data->bcf_header);
	hts_idx_destroy(data->bcf_index);
	bcf_close(data->bcf_file);
}

void BCFReaderPop(WiggleIterator * wi) {
	BCFReaderData * data = (BCFReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void BcfReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BCFReaderData * data = (BCFReaderData *) wi->data;

	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}
	// Note: BCF encoding is 0 based, hence -1s
	data->bcf_iterator = bcf_itr_queryi(data->bcf_index, bcf_hdr_name2id(data->bcf_header, chrom), start - 1, finish - 1);
	if (data->bcf_iterator == NULL) {
		fprintf(stderr, "Could not find index file to BCF file %s.\n", data->filename);
		exit(1);
	}
	launchBufferedReader(&downloadBCFFile, data, &(data->bufferedReaderData));
	wi->done = false;
	BCFReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BCFReaderPop(wi);
}

WiggleIterator * BcfReader(char * filename, bool holdFire) {
	BCFReaderData * data = (BCFReaderData *) calloc(1, sizeof(BCFReaderData));
	OpenBCFFile(data, filename);
	if (!holdFire)
		launchBufferedReader(&downloadBCFFile, data, &(data->bufferedReaderData));
	return newWiggleIterator(data, &BCFReaderPop, &BcfReaderSeek, 0, true);
}
