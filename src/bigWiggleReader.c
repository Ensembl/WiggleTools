// Copyright 2013-2014 EMBL-EBI
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

bool readBigWigBuffer(BigFileReaderData * data) {
	int i;
	int start, finish;
	double value;
	char *blockPt = data->uncompressBuf;
	struct bwgSectionHead head;
	start = finish = 0;
	
	bwgSectionHeadFromMem(&(blockPt), &(head), data->isSwapped);

        for (i = 0; i < head.itemCount; i++) {
		switch (head.type)
		    {
		    case bwgTypeBedGraph:
			{
			// +1 because BigWig coords are 0-based...
			start = memReadBits32(&(blockPt), data->isSwapped) + 1;
			finish = memReadBits32(&(blockPt), data->isSwapped) + 1;
			value = memReadFloat(&(blockPt), data->isSwapped);
			break;
			}
		    case bwgTypeVariableStep:
			{
			// +1 because BigWig coords are 0-based...
			start = memReadBits32(&(blockPt), data->isSwapped) + 1;
			finish = start + head.itemSpan;
			value = memReadFloat(&(blockPt), data->isSwapped);
			break;
			}
		    case bwgTypeFixedStep:
			{
			if (i==0) 
			    {
			    // +1 because BigWig coords are 0-based...
			    start = head.start + 1;
			    finish = start + head.itemSpan;
			    value = memReadFloat(&(blockPt), data->isSwapped);
			    }
			else
			    {
			    start += head.itemStep;
			    finish += head.itemStep;
			    value = memReadFloat(&(blockPt), data->isSwapped);
			    }
			break;
			}
		    default:
			{
			fprintf(stderr, "Unrecognized head type in Wiggle file\n");
			exit(1);
			}
		    }
		if (data->stop > 0) {
			if (start >= data->stop) {
				return true;
			} else if (finish > data->stop) {
				finish = data->stop;
			}
		}

		if (pushValuesToBuffer(data->bufferedReaderData, data->chrom, start, finish, value))
			return true;
	}

	return false;
}

static void openBigWigFile(BigFileReaderData * data, char * filename, bool holdFire) {
	data->filename = filename;
	data->bwf = bigWigFileOpen(data->filename);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;
	data->readBuffer = &readBigWigBuffer;
	data->uncompressBuf = (char *) needLargeMem(data->bwf->uncompressBufSize);
	if (!holdFire)
		launchBufferedReader(&downloadBigFile, data, &(data->bufferedReaderData));
}

WiggleIterator * BigWiggleReader(char * f, bool holdFire) {
	BigFileReaderData * data = (BigFileReaderData *) calloc(1, sizeof(BigFileReaderData));
	openBigWigFile(data, f, holdFire);
	return newWiggleIterator(data, &BigFileReaderPop, &BigFileReaderSeek, 0);
}	
