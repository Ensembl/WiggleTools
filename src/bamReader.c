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
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "wiggleIterator.h"
#include "bufferedReader.h"
#include "hashfib.h"

typedef struct bamFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int chrom_tid;
	int start, stop;
	bool read_count;
	BufferedReaderData * bufferedReaderData;

	// BAM stuff
	samFile * fp;
	hts_idx_t * idx;
	bam_hdr_t * header;
} BamReaderData;

static void storeReadComponents(HashFib * starts, HashFib * ends, bam1_t * aln) {
	// Sometimes a read has coordinates, but is not mapped (cf BWA)
	// We should skip these exceptions
	if (aln->core.flag & 0x4)
		return;
	// Note that BAM coords are 0-based, hence +1
	int start = aln->core.pos + 1;
	uint32_t *cigar = bam_get_cigar(aln);
	int k;
	for (k = 0; k < aln->core.n_cigar; ++k) {
		int operator = cigar[k]&BAM_CIGAR_MASK;
		int length = cigar[k]>>BAM_CIGAR_SHIFT;
		switch (operator) {
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
			case BAM_CDEL:
				hashfib_insert(starts, start);
				hashfib_insert(ends, start + length);
			case BAM_CREF_SKIP:
				start += length;
		}
	}
}

static bam1_t * nextRead(BamReaderData * data, hts_itr_t * iter, bam1_t * aln) {
	if (sam_itr_next(data->fp, iter, aln) < 0) {
		bam_destroy1(aln);
		return NULL;
	} else  {
		return aln;
	}
}

static bool consumeIteratorWithCigars(BamReaderData * data, hts_itr_t * iter, char * query_chrom, int query_chrom_tid, int query_stop) {
	// Iterate through data
	bam1_t *aln = bam_init1();
	HashFib * starts = hashfib_construct();
	HashFib * ends = hashfib_construct();
	int start, finish = -1, chrom_tid = -1;
	int value = 0;

	aln = nextRead(data, iter, aln);

	while(1) {
		// Remove dead weight
		if (!hashfib_empty(ends) && hashfib_min(ends) == finish)
			value -= hashfib_remove_min(ends);

		// Stream more data into heaps
		while (aln
		       && (aln->core.tid == chrom_tid || hashfib_empty(ends)) 
		       && (hashfib_empty(ends) || aln->core.pos <= hashfib_min(ends) || hashfib_empty(starts) || aln->core.pos <= hashfib_min(starts))
		) {
			storeReadComponents(starts, ends, aln);
			chrom_tid = aln->core.tid;
			aln = nextRead(data, iter, aln);
		}

		if (!hashfib_empty(ends)) {
			// Choose new start
			if (value)
				start = finish;
			else
				start = hashfib_min(starts);

			// If gone overboard
			if ((start >= query_stop && chrom_tid == query_chrom_tid) || chrom_tid > query_chrom_tid) {
				hashfib_destroy(starts);
				hashfib_destroy(ends);
				return false;
			}

			// Load new weight
			if (!hashfib_empty(starts) && hashfib_min(starts) == start) 
				value += hashfib_remove_min(starts);

			// Choose finish
			if (hashfib_empty(starts) || hashfib_min(ends) < hashfib_min(starts))
				finish = hashfib_min(ends);
			else 
				finish = hashfib_min(starts);

			// If gone overboard
			if (finish > query_stop)
				finish = query_stop;

			// Push out
			if (pushValuesToBuffer(data->bufferedReaderData, data->header->target_name[chrom_tid], start, finish, value)) {
				hashfib_destroy(starts);
				hashfib_destroy(ends);
				return true;
			}
		} else {
			// No ends => end of file iterator
			hashfib_destroy(starts);
			hashfib_destroy(ends);
			return false;
		}
	}
}

static bool consumeIteratorNoCigars(BamReaderData * data, hts_itr_t * iter, char * query_chrom, int query_chrom_tid, int query_stop) {
	// Iterate through data
	bam1_t *aln = bam_init1();
	int start, finish, value;

	aln = nextRead(data, iter, aln);

	while(aln && aln->core.tid == query_chrom_tid && aln->core.pos < query_stop) {
		// Read coordinates
		// Note that BAM coords are 0-based, hence +1
		start = aln->core.pos + 1;
		finish = start + 1;

		// Compute value
		value = 0;
		while (aln
		       && aln->core.tid == query_chrom_tid 
		       && aln->core.pos == start - 1) {
			value++;
			aln = nextRead(data, iter, aln);
		}

		// Push on
		if (pushValuesToBuffer(data->bufferedReaderData, query_chrom, start, finish, value)) 
			return true;
	}

	return false;
}

static bool downloadBamFileChromosome(BamReaderData * data, char * query_chrom, uint32_t query_start, uint32_t query_stop) {
	// Create iterator
	int query_chrom_tid = bam_name2id(data->header, query_chrom);
	hts_itr_t * iter = sam_itr_queryi(data->idx, query_chrom_tid, query_start, query_stop);
	if(data->header == NULL || iter == NULL) {
		fprintf(stderr, "Unable to iterate to region %s:%i-%i within %s BAM file.", query_chrom, query_start, query_stop, data->filename);
		exit(1);
	}
	bool status;
	if (data->read_count) 
		status = consumeIteratorNoCigars(data, iter, query_chrom, query_chrom_tid, query_stop);
	else 
		status = consumeIteratorWithCigars(data, iter, query_chrom, query_chrom_tid, query_stop);

	if (iter)
		bam_itr_destroy(iter);
	return status;
}

typedef struct {
	char * name;
	int tid;
} name_id_st;

static int comp_name_id_st(const void * A, const void * B) {
	char * A_name = ((name_id_st *) A)->name;
	char * B_name = ((name_id_st *) B)->name;
	return strcmp(A_name, B_name);
}

static void * downloadBamFile(void * args) {
	BamReaderData * data = (BamReaderData *) args;

	if (data->chrom) {
		downloadBamFileChromosome(data, data->chrom, data->start, data->stop);
	} else {
		// Copy chromosome name into array of name_id_st:
		name_id_st * name_ids = calloc(data->header->n_targets, sizeof(name_id_st));
		int index;
		for (index = 0; index < data->header->n_targets; index++) {
			name_ids[index].name = data->header->target_name[index];
			name_ids[index].tid = index;
		}

		// Sort list of chromosomes alphabetically:
		qsort(name_ids, data->header->n_targets, sizeof(name_id_st), comp_name_id_st); 
		
		// Iterate through chromsomes in alphabetic order
		for (index = 0; index < data->header->n_targets; index++)
			if (downloadBamFileChromosome(data, name_ids[index].name, 0, data->header->target_len[name_ids[index].tid]))
				break;
	}

	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void OpenBamFile(BamReaderData * data, char * filename) {
	data->filename = filename;

	// read the header and initialize data
	data->fp = hts_open(filename, "r");

        //Get the header
        data->header = sam_hdr_read(data->fp);

	// Get the index
	data->idx = sam_index_load(data->fp, data->filename);
	if(data->idx == NULL) {
		fprintf(stderr, "Unable to open BAM/SAM index. Make sure alignments are indexed\n");
		exit(1);
	}
}

void closeBamFile(BamReaderData * data) {
	if (data->idx)
		hts_idx_destroy(data->idx);
        bam_hdr_destroy(data->header);
        sam_close(data->fp);
}

void BamReaderPop(WiggleIterator * wi) {
	BamReaderData * data = (BamReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void BamReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BamReaderData * data = (BamReaderData *) wi->data;

	// Kill ongoing jobs
	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}

	// Set boundaries
	data->chrom = chrom;
	data->start = start;
	data->stop = finish;

	// Weeeee
	launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	wi->done = false;
	BamReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(wi->chrom, chrom) == 0 && wi->finish <= start)))
		BamReaderPop(wi);

}

WiggleIterator * BamReader(char * filename, bool holdFire, bool read_count) {
	BamReaderData * data = (BamReaderData *) calloc(1, sizeof(BamReaderData));
	data->read_count = read_count;
	OpenBamFile(data, filename);
	if (!holdFire)
		launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek, 0, false);
}
