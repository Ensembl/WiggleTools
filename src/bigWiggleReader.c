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
		killDownloader(data);
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
	if (data->stop > 0) {
		if (wi->start >= data->stop) {
			killDownloader(data);
			wi->done = true;
			return;
		} else if (wi->finish > data->stop) {
			wi->finish = data->stop;
		}
	}

        if (++(data->i) == data->head.itemCount)
	    BigWiggleReaderGoToNextBlock(data);

}

void BigWiggleReaderSeek(WiggleIterator * wi, const char * chrom, int start, int finish) {
	BigFileReaderData * data = (BigFileReaderData *) wi->data; 

	killDownloader(data);
	data->chrom = chrom;
	data->stop = finish;
	wi->done = false;
	launchDownloader(data);
	BigWiggleReaderEnterBlock(data);
	BigWiggleReaderPop(wi);

	while (!wi->done && (strcmp(wi->chrom, chrom) < 0 || (strcmp(chrom, wi->chrom) == 0 && wi->finish <= start)))
		BigWiggleReaderPop(wi);

	if (!wi->done && strcmp(chrom, wi->chrom) == 0 && wi->start < start)
		wi->start = start;
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
	BigWiggleReaderEnterBlock(data);

	return newWiggleIterator(data, &BigWiggleReaderPop, &BigWiggleReaderSeek);
}	
