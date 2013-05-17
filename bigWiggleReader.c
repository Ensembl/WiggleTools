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
#include "wiggleTools.h"
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
	/* File variables */
	struct bbiFile* bwf;
	boolean isSwapped;

	// Internal loop variables:
	// Chromosomes
	struct bbiChromInfo *chrom, *chromList;

	// Sets of blocks within chromosome
	struct fileOffsetSize *blockList, *block, *beforeGap, *afterGap;
	struct udcFile *udc;
	char *uncompressBuf;

	// Blocks within set of blocks
        bits64 mergedOffset, mergedSize;
        char *mergedBuf, *blockBuf;
	char *blockPt, *blockEnd;
	struct bwgSectionHead head;

	// Items within block
	bits16 i;
} BigWiggleReaderData;

void BigWiggleReaderEnterBlock(BigWiggleReaderData * data) {
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

void BigWiggleReaderEnterSetOfBlocks(BigWiggleReaderData * data) {
       /* Find contigious blocks and read them into mergedBuf. */
       fileOffsetSizeFindGap(data->block, &(data->beforeGap), &(data->afterGap));
       data->mergedOffset = data->block->offset;
       data->mergedSize = data->beforeGap->offset + data->beforeGap->size - data->mergedOffset;
       udcSeek(data->udc, data->mergedOffset);
       data->mergedBuf = (char *) needLargeMem(data->mergedSize);
       udcMustRead(data->udc, data->mergedBuf, data->mergedSize);
       data->blockBuf = data->mergedBuf;
       BigWiggleReaderEnterBlock(data); 
}

void BigWiggleReaderEnterChromosome(BigWiggleReaderData * data) {
	data->block = data->blockList = bbiOverlappingBlocks(data->bwf, data->bwf->unzoomedCir, data->chrom->name, 0, data->chrom->size, NULL);
	BigWiggleReaderEnterSetOfBlocks(data);
}

void BigWiggleReaderOpenFile(BigWiggleReaderData * data, char * f) {
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

void BigWiggleReaderGoToNextChromosome(BigWiggleReaderData * data) {
	data->chrom = data->chrom->next;
	if (!data->chrom) {
	    bbiFileClose(&(data->bwf));
	} else 
	    BigWiggleReaderEnterChromosome(data);
}

void BigWiggleReaderGoToNextSetOfBlocks(BigWiggleReaderData * data) {
	if (data->block == NULL) {
		freeMem(data->uncompressBuf);
		BigWiggleReaderGoToNextChromosome(data);
	} else 
		BigWiggleReaderEnterSetOfBlocks(data);
}

void BigWiggleReaderGoToNextBlock(BigWiggleReaderData * data) {
	data->blockBuf += data->block->size;
	data->block = data->block->next;
	if (data->block == data->afterGap) {
		freeMem(data->mergedBuf);
		BigWiggleReaderGoToNextSetOfBlocks(data);
	} else 
		BigWiggleReaderEnterBlock(data);
}

void BigWiggleReaderPop(WiggleIterator * wi) {
	BigWiggleReaderData * data;

	if (wi->done)
		return;

	data = (BigWiggleReaderData*) wi->data;

	if (!data->chrom) {
		// Because strings are passed by reference instead of pointers deallocating
		// These objects would mess up objects for iterators downstream
		//bbiChromInfoFreeList(&(data->chromList));
		wi->done = true;
		return;
	} else {
		wi->chrom = data->chrom->name;
	}

	switch (data->head.type)
	    {
	    case bwgTypeBedGraph:
		{
		wi->start = memReadBits32(&(data->blockPt), data->isSwapped);
		wi->finish = memReadBits32(&(data->blockPt), data->isSwapped);
		wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		break;
		}
	    case bwgTypeVariableStep:
		{
		wi->start = memReadBits32(&(data->blockPt), data->isSwapped);
		wi->finish = wi->start + data->head.itemSpan;
		wi->value = memReadFloat(&(data->blockPt), data->isSwapped);
		break;
		}
	    case bwgTypeFixedStep:
		{
		if (data->i==0) 
		    {
		    wi->start = data->head.start;
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
