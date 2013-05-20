// Copyright (c) 2013, Daniel Zerbino
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
// (1) Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer. 
// 
// (2) Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in
// the documentation and/or other materials provided with the
// distribution.  
// 
// (3)The name of the author may not be used to
// endorse or promote products derived from this software without
// specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "stdio.h"

// Local header
#include "wiggleIterators.h"

// Kent library headers
#include "common.h"
#include "zlibFace.h"
#include "bbiFile.h"
#include "bigWig.h"
#include "udc.h"
#include "bwgInternal.h"

//////////////////////////////////////////////////////
// Compressed File Reader
//////////////////////////////////////////////////////

typedef struct bigWiggleReaderData_st {
	// File level variables
	struct bbiFile* bwf;
	struct udcFile *udc;
	boolean isSwapped;
	char *uncompressBuf;
	struct bbiChromInfo *chromList;

	// Chromosome within file
	struct bbiChromInfo *chrom;
	struct fileOffsetSize *blockList;

	// Contiguous runs of blocks within chromosome
	struct fileOffsetSize *afterGap;
	char * mergedBuf;

	// Blocks within run of blocks
	struct fileOffsetSize * block;
        char *blockBuf, *blockEnd;
	struct bwgSectionHead head;

	// Items within block
	char *blockPt;
	bits16 i;
} BigWiggleReaderData;

static void BigWiggleReaderEnterBlock(BigWiggleReaderData * data) {
	/* Uncompress if necessary. */
	if (data->uncompressBuf)
	    {
	    data->blockPt = data->uncompressBuf;
	    int uncSize = zUncompress(data->blockBuf, data->block->size, data->uncompressBuf, data->bwf->uncompressBufSize);
	    data->blockEnd = data->blockPt + uncSize;
	    }
	else
	    {
	    data->blockPt = data->blockBuf;
	    data->blockEnd = data->blockPt + data->block->size;
	    }

	/* Deal with insides of block. */
	bwgSectionHeadFromMem(&(data->blockPt), &(data->head), data->isSwapped);

	data->i = 0;
}

static void BigWiggleReaderEnterRunOfBlocks(BigWiggleReaderData * data) {
       /* Find contigious blocks and read them into mergedBuf. */
       struct fileOffsetSize *beforeGap;
       bits64 mergedSize;
       fileOffsetSizeFindGap(data->block, &beforeGap, &(data->afterGap));
       mergedSize = beforeGap->offset + beforeGap->size - data->block->offset;
       udcSeek(data->udc, data->block->offset);
       data->mergedBuf = (char *) needLargeMem(mergedSize);
       udcMustRead(data->udc, data->mergedBuf, mergedSize);
       data->blockBuf = data->mergedBuf;
       BigWiggleReaderEnterBlock(data); 
}

static void BigWiggleReaderEnterChromosome(BigWiggleReaderData * data) {
	data->block = data->blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, data->chrom->name, 0, data->chrom->size, NULL);
	BigWiggleReaderEnterRunOfBlocks(data);
}

static void BigWiggleReaderOpenFile(BigWiggleReaderData * data, char * f) {
	/* Opening up BigWig file */
	data->bwf = bigWigFileOpen(f);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;

	/* Set up for uncompression optionally. */
	if (data->bwf->uncompressBufSize > 0)
		data->uncompressBuf = (char *) needLargeMem(data->bwf->uncompressBufSize);
	else
		data->uncompressBuf = NULL;

	data->chrom = data->chromList = bbiChromList(data->bwf);
	BigWiggleReaderEnterChromosome(data);
}

static void BigWiggleReaderCloseFile(BigWiggleReaderData * data) {
	// Because strings are passed by reference instead of pointers deallocating
	// this linked list would mess up objects for iterators downstream
	//bbiChromInfoFreeList(&(data->chromList));
	freeMem(data->uncompressBuf);
	bbiFileClose(&(data->bwf));
}

static void BigWiggleReaderGoToNextChromosome(BigWiggleReaderData * data) {
	slFreeList(&(data->blockList));
	if ((data->chrom = data->chrom->next))
		BigWiggleReaderEnterChromosome(data);
	else 
		BigWiggleReaderCloseFile(data);
}

static void BigWiggleReaderGoToNextRunOfBlocks(BigWiggleReaderData * data) {
	freeMem(data->mergedBuf);
	if (data->block)
		BigWiggleReaderEnterRunOfBlocks(data);
	else 
		BigWiggleReaderGoToNextChromosome(data);
}

static void BigWiggleReaderGoToNextBlock(BigWiggleReaderData * data) {
	data->blockBuf += data->block->size;
	if ((data->block = data->block->next) == data->afterGap)
		BigWiggleReaderGoToNextRunOfBlocks(data);
	else 
		BigWiggleReaderEnterBlock(data);
}

void BigWiggleReaderPop(WiggleIterator * wi) {
	BigWiggleReaderData * data;

	if (wi->done)
		return;

	data = (BigWiggleReaderData*) wi->data;

	if (!data->chrom) {
		// Passive agressive indicator that the iterator should be closed
		// (avoids passing the WiggleIterator reference needlessly across 
		// all functions).
		wi->done = true;
		return;
	} else {
		wi->chrom = data->chrom->name;
	}

	switch (data->head.type)
	    {
	    case bwgTypeBedGraph:
		{
		// +1 because BigWig coords are 0-based...
		wi->start = memReadBits32(&(data->blockPt), data->isSwapped) + 1;
		wi->finish = memReadBits32(&(data->blockPt), data->isSwapped) + 1;
		wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		break;
		}
	    case bwgTypeVariableStep:
		{
		// +1 because BigWig coords are 0-based...
		wi->start = memReadBits32(&(data->blockPt), data->isSwapped) + 1;
		wi->finish = wi->start + data->head.itemSpan;
		wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		break;
		}
	    case bwgTypeFixedStep:
		{
		if (data->i==0) 
		    {
	 	    // +1 because BigWig coords are 0-based...
		    wi->start = data->head.start + 1;
		    wi->finish = wi->start + data->head.itemSpan;
		    wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		    }
		else
		    {
		    wi->start += data->head.itemStep;
		    wi->finish += data->head.itemStep;
		    wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		    }
		break;
		}
	    default:
		{
		fprintf(stderr, "Unrecognized head type in Wiggle file\n");
		exit(1);
		break;
		}
	    }
        if (++(data->i) == data->head.itemCount)
	    {
	    assert(data->blockPt == data->blockEnd);
	    BigWiggleReaderGoToNextBlock(data);
	    }

}

WiggleIterator * BigWiggleReader(char * f) {
	BigWiggleReaderData * data = (BigWiggleReaderData *) calloc(1, sizeof(BigWiggleReaderData));
	BigWiggleReaderOpenFile(data, f);
	WiggleIterator * new = newWiggleIterator(data, &BigWiggleReaderPop);
	return new;
}	
