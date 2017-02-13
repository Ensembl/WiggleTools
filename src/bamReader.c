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
#include "fib.h"

typedef struct bamFileReaderData_st {
	// Arguments to downloader
	char * filename;
	char * chrom;
	int chrom_tid;
	int start, stop;
	BufferedReaderData * bufferedReaderData;

	// BAM stuff
	samFile * fp;
	hts_idx_t * idx;
	bam_hdr_t * header;
} BamReaderData;

static void storeReadComponents(FibHeap * starts, FibHeap * ends, bam1_t * aln) {
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
				fh_insert(starts, start, 0);
				fh_insert(ends, start + length, 0);
			case BAM_CREF_SKIP:
				start += length;
		}
	}
}

static bam1_t * nextRead(BamReaderData * data, hts_itr_t * iter, bam1_t * aln) {
	if ((iter && sam_itr_next(data->fp, iter, aln) < 0) 
	     || (!iter && sam_read1(data->fp, data->header, aln) < 0)
	   ) {
		bam_destroy1(aln);
		return NULL;
	} else  {
		return aln;
	}
}

static void * downloadBamFile(void * args) {
	BamReaderData * data = (BamReaderData *) args;

	//Initialize iterator
        hts_itr_t *iter = NULL;
	if (data->chrom) {
		data->chrom_tid = bam_name2id(data->header, data->chrom);
		iter = sam_itr_queryi(data->idx, data->chrom_tid, data->start, data->stop);
		if(data->header == NULL || iter == NULL) {
			fprintf(stderr, "Unable to iterate to region within BAM.");
			exit(1);
		}
	}

	// Iterate through data
	bam1_t *aln = bam_init1();
	FibHeap * starts = fh_makeheap();
	FibHeap * ends = fh_makeheap();
	int start, finish = -1, chrom_tid = -1;
	int value = 0;

	aln = nextRead(data, iter, aln);

	while(1) {
		// Remove dead weight
		while (fh_notempty(ends) && fh_min(ends) == finish) {
			value--;
			fh_extractmin(ends);
		}

		// Stream more data into heaps
		while (aln
		       && (aln->core.tid == chrom_tid || fh_empty(ends)) 
		       && (fh_empty(ends) || aln->core.pos <= fh_min(ends) || fh_empty(starts) || aln->core.pos <= fh_min(starts))
		) {
			storeReadComponents(starts, ends, aln);
			chrom_tid = aln->core.tid;
			aln = nextRead(data, iter, aln);
		}

		if (fh_notempty(ends)) {
			// Choose new start
			if (value)
				start = finish;
			else
				start = fh_min(starts);

			// Load new weight
			while (fh_notempty(starts) && fh_min(starts) == start) {
				value++;
				fh_extractmin(starts);
			}

			// Choose finish
			if (fh_notempty(starts))
				finish = fh_min(starts);
			if (fh_notempty(ends) && (fh_empty(starts) || fh_min(ends) < finish))
				finish = fh_min(ends);

			// Possible cutoff of data
			if (data->stop > 0) {
				if ((start >= data->stop && chrom_tid == data->chrom_tid) || chrom_tid > data->chrom_tid) {
					fh_deleteheap(starts);
					fh_deleteheap(ends);
					if (iter)
						bam_itr_destroy(iter);
					endBufferedSignal(data->bufferedReaderData);
					return NULL;
				} else if (finish > data->stop)
					finish = data->stop;
			}

			if (pushValuesToBuffer(data->bufferedReaderData, data->header->target_name[chrom_tid], start, finish, value)) {
				fh_deleteheap(starts);
				fh_deleteheap(ends);
				if (iter)
					bam_itr_destroy(iter);
				return NULL;
			}
		} else {
			// No ends => end of file iterator
			fh_deleteheap(starts);
			fh_deleteheap(ends);
			if (iter)
				bam_itr_destroy(iter);
			endBufferedSignal(data->bufferedReaderData);
			return NULL;
		}
	}
}

void OpenBamFile(BamReaderData * data, char * filename) {
	data->filename = filename;

	// read the header and initialize data
	data->fp = hts_open(filename, "r");

        //Get the header
        data->header = sam_hdr_read(data->fp);
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

	//Load the index
	if (!data->idx) {
		data->idx = sam_index_load(data->fp, data->filename);
		if(data->idx == NULL) {
			fprintf(stderr, "Unable to open BAM/SAM index. Make sure alignments are indexed\n");
			exit(1);
		}
	}

	// Kill ongoing jobs
	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}

	// Set boundaries
	data->chrom = chrom;
	data->stop = finish;

	// Weeeee
	launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	wi->done = false;
	BamReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(wi->chrom, chrom) == 0 && wi->finish <= start)))
		BamReaderPop(wi);

}

WiggleIterator * BamReader(char * filename, bool holdFire) {
	BamReaderData * data = (BamReaderData *) calloc(1, sizeof(BamReaderData));
	OpenBamFile(data, filename);
	if (!holdFire)
		launchBufferedReader(&downloadBamFile, data, &(data->bufferedReaderData));
	return newWiggleIterator(data, &BamReaderPop, &BamReaderSeek, 0);
}
