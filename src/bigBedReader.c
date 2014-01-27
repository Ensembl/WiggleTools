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

#include "bigFileReader.h"

static void BigBedReaderEnterBlock(BigFileReaderData * data) {
	if (data->blockData)
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
		killDownloader(data);
		wi->done = true;
		return;
	}

	/* Read next record into local variables. */
	memReadBits32(&data->blockPt, data->isSwapped);	// Read and discard chromId
	wi->chrom = data->chrom;
	wi->start = memReadBits32(&data->blockPt, data->isSwapped) + 1;
	wi->finish = memReadBits32(&data->blockPt, data->isSwapped) + 1; 

	// Skip boring stuff...
	for (;;)
		if (*(data->blockPt++) <= 0)
			break;

	if (data->stop > 0) {
		if (wi->start >= data->stop) {
			killDownloader(data);
			wi->done = true;
			return;
		} else if (wi->finish > data->stop) {
			wi->finish = data->stop;
		}
	}

	if (data->blockPt == data->blockData->blockEnd)
		BigBedReaderGoToNextBlock(data);
}

void BigBedReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data; 

	killDownloader(data);
	data->chrom = chrom;
	data->stop = finish;
	wi->done = false;
	launchDownloader(data);
	BigBedReaderEnterBlock(data);
	BigBedReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start))) 
		BigBedReaderPop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
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
	WiggleIterator * res = newWiggleIterator(data, &BigBedReaderPop, &BigBedReaderSeek);
	res->overlaps = true;
	return res;
}	
