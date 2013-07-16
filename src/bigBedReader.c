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

#include "bigFileReader.h"

static void BigBedReaderEnterBlock(BigFileReaderData * data) {
	if (!data->blockData)
		abort();
	enterBlock(data);
}

void BigBedReaderGoToNextBlock(BigFileReaderData * data) {
	goToNextBlock(data);
	if (data->blockData)
		BigBedReaderEnterBlock(data);
}

void BigBedReaderPop(WiggleIterator * wi) {
	if (wi->done)
		return;

	BigFileReaderData * data = (BigFileReaderData*) wi->data;

	if (!data->blockData) {
		wi->done = true;
		return;
	}

	/* Read next record into local variables. */
	memReadBits32(&data->blockPt, data->isSwapped);	// Read and discard chromId
	wi->chrom = data->chrom;
	wi->start = memReadBits32(&data->blockPt, data->isSwapped);
	wi->finish = memReadBits32(&data->blockPt, data->isSwapped) + 1; 

	// Skip boring stuff...
	for (;;)
		if (*(data->blockPt++) <= 0)
			break;

	if (data->stop > 0 && wi->start > data->stop) {
		wi->done = true;
		return;
	}

	if (data->blockPt == data->blockData->blockEnd)
		BigBedReaderGoToNextBlock(data);
}

void BigBedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data; 

	killDownloader(data);
	data->chrom = chrom;
	data->stop = finish;
	launchDownloader(data);
	enterBlock(data);
	wi->done = false;
	BigBedReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigBedReaderPop(wi);
}

static void openBigBedFile(BigFileReaderData * data) {
	data->bwf = bigBedFileOpen(data->filename);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;
}

WiggleIterator * BigBedReader(char * f) {
	BigFileReaderData * data = (BigFileReaderData *) calloc(1, sizeof(BigFileReaderData));
	data->filename = f;
	openBigBedFile(data);
	launchDownloader(data);
	BigBedReaderEnterBlock(data);
	return UnionWiggleIterator(newWiggleIterator(data, &BigBedReaderPop, &BigBedReaderSeek));
}	
