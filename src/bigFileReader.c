// Copyright [1999-2016] EMBL-European Bioinformatics Institute
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

#include "bigFileReader.h"
#include "bufferedReader.h"

static int MAX_BLOCKS = 100;

static bool openBlock(BigFileReaderData * data, struct fileOffsetSize * block, char * blockBuf) {
	size_t uncompressBufSize = data->bwf->uncompressBufSize;

	if (uncompressBufSize > 0) {
		int uncSize = zUncompress(blockBuf, block->size, data->uncompressBuf, uncompressBufSize);
		data->blockEnd = data->uncompressBuf + uncSize;
	} else {
		data->uncompressBuf = blockBuf;
		data->blockEnd = blockBuf + block->size;
	}

	// Callback to specialised BigBed or BigWig function
	return data->readBuffer(data);
}

static bool downloadBlockRun(BigFileReaderData * data, char * chrom, struct fileOffsetSize * firstBlock, struct fileOffsetSize * afterBlock, bits64 mergedSize) {
	char * mergedBuf, *blockBuf;
	struct fileOffsetSize * block;

	udcSeek(data->udc, firstBlock->offset);
	blockBuf = mergedBuf = (char *) needLargeMem(mergedSize);
	udcMustRead(data->udc, mergedBuf, mergedSize);

	for (block = firstBlock; block != afterBlock; block = block->next) {
		if (openBlock(data, block, blockBuf)) {
			freeMem(mergedBuf);
			return true;
		}
		blockBuf += block->size;
	}

	freeMem(mergedBuf);
	return false;
}

static bool downloadBigRegion(BigFileReaderData * data, char * chrom, int start, int finish) {
	struct fileOffsetSize *blockList, *block, *beforeGap, *afterGap;
	int blockCounter;
	bits64 mergedSize;

	data->chrom = chrom;
	blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, chrom, start, finish, NULL);

	for (block = blockList; block; block=afterGap) {
		/* Read contiguous blocks into mergedBuf. */
		fileOffsetSizeFindGap(block, &beforeGap, &afterGap);

		// Little hack to limit the number of blocks read at any time
		struct fileOffsetSize * blockPtr, * prevBlock;
		blockCounter = 0;
		prevBlock = block;
		for (blockPtr = block; blockPtr != afterGap && blockCounter < MAX_BLOCKS; blockPtr = blockPtr->next) {
			blockCounter++;
			prevBlock = blockPtr;
		}
		if (blockCounter == MAX_BLOCKS) {
			beforeGap = prevBlock;
			afterGap = blockPtr;
		}

		mergedSize = beforeGap->offset + beforeGap->size - block->offset;

		if (downloadBlockRun(data, chrom, block, afterGap, mergedSize)) {
			slFreeList(blockList);
			return true;
		}
	}

	if (blockList)
		slFreeList(blockList);
	return false;
}

static void downloadFullGenome(BigFileReaderData * data) {
	struct bbiChromInfo *chromList = bbiChromList(data->bwf);
	struct bbiChromInfo *chrom;

	for (chrom = chromList; chrom; chrom = chrom->next) 
		if (downloadBigRegion(data, chrom->name, 0, chrom->size))
			break;

	// TODO free chromList memory... yes but labels lost! Need to be copied out first....
	//bbiChromInfoFreeList(&(data->chromList));
}

void * downloadBigFile(void * args) {
	BigFileReaderData * data = (BigFileReaderData *) args;

	if (!data->chrom)
		downloadFullGenome(data);
	else 
		downloadBigRegion(data, data->chrom, data->start, data->stop);

	endBufferedSignal(data->bufferedReaderData);
	return NULL;
}

void BigFileReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data; 

	if (data->bufferedReaderData) {
		killBufferedReader(data->bufferedReaderData);
		free(data->bufferedReaderData);
		data->bufferedReaderData = NULL;
	}
	data->chrom = chrom;
	data->start = start;
	data->stop = finish;
	launchBufferedReader(&downloadBigFile, data, &(data->bufferedReaderData));
	wi->done = false;
	BigFileReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigFileReaderPop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;

}

void BigFileReaderPop(WiggleIterator * wi) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data;
	BufferedReaderPop(wi, data->bufferedReaderData);
}

void BigFileReaderCloseFile(BigFileReaderData * data) {
	bbiFileClose(&(data->bwf));
}
