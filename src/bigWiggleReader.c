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

// Local header
#include "bigFileReader.h"

static void BigWiggleReaderEnterBlock(BigFileReaderData * data) {
	if (data->blockData) {
		enterBlock(data);
		bwgSectionHeadFromMem(&(data->blockPt), &(data->head), data->isSwapped);
		data->i = 0;
	}
}

void BigWiggleReaderGoToNextBlock(BigFileReaderData * data) {
	goToNextBlock(data);
	BigWiggleReaderEnterBlock(data);
}

void BigWiggleReaderPop(WiggleIterator * wi) {
	BigFileReaderData * data;

	if (wi->done)
		return;

	data = (BigFileReaderData*) wi->data;

	if (!data->blockData) {
		wi->done = true;
		return;
	}
		
	wi->chrom = data->chrom;

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
		killDownloader(data);
		exit(1);
		break;
		}
	    }
	if (data->stop > 0 && (wi->start > data->stop)) {
		wi->done = true;
		return;
	}

        if (++(data->i) == data->head.itemCount)
	    BigWiggleReaderGoToNextBlock(data);

}

void BigWiggleReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data; 

	killDownloader(data);
	data->chrom = chrom;
	data->stop = finish;
	launchDownloader(data);
	if (data->blockData)
		BigWiggleReaderEnterBlock(data);
	wi->done = false;
	BigWiggleReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigWiggleReaderPop(wi);
}

static void openBigWigFile(BigFileReaderData * data) {
	data->bwf = bigWigFileOpen(data->filename);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;
}

WiggleIterator * BigWiggleReader(char * f) {
	BigFileReaderData * data = (BigFileReaderData *) calloc(1, sizeof(BigFileReaderData));
	data->filename = f;
	openBigWigFile(data);
	launchDownloader(data);
	if (data->blockData)
		BigWiggleReaderEnterBlock(data);

	return newWiggleIterator(data, &BigWiggleReaderPop, &BigWiggleReaderSeek);
}	
