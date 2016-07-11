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

bool readBigBedBuffer(BigFileReaderData * data) {
	char *blockPt;

	/* Read next record into local variables. */
	for (blockPt = data->uncompressBuf; blockPt != data->blockEnd; ) {
		memReadBits32(&blockPt, data->isSwapped);	// Read and discard chromId
		int start = memReadBits32(&blockPt, data->isSwapped) + 1;
		int finish = memReadBits32(&blockPt, data->isSwapped) + 1; 

		// Skip boring stuff...
		for (;;)
			if (*(blockPt++) <= 0)
				break;

		if (data->stop > 0) {
			if (start >= data->stop)
				return true;
			else if (finish > data->stop)
				finish = data->stop;
		}
		if (pushValuesToBuffer(data->bufferedReaderData, data->chrom, start, finish, 1))
			return true;
	}

	return false;
}

void openBigBedFile(BigFileReaderData * data, char * filename, bool holdFire) {
	data->filename = filename;
	data->bwf = bigBedFileOpen(data->filename);
	data->isSwapped = data->bwf->isSwapped;
	bbiAttachUnzoomedCir(data->bwf);
	data->udc = data->bwf->udc;
	data->readBuffer = &readBigBedBuffer;
	data->uncompressBuf = (char *) needLargeMem(data->bwf->uncompressBufSize);
	if (!holdFire)
		launchBufferedReader(&downloadBigFile, data, &(data->bufferedReaderData));
}

WiggleIterator * BigBedReader(char * f, bool holdFire) {
	BigFileReaderData * data = (BigFileReaderData *) calloc(1, sizeof(BigFileReaderData));
	openBigBedFile(data, f, holdFire);
	WiggleIterator * res = newWiggleIterator(data, &BigFileReaderPop, &BigFileReaderSeek, 0);
	res->overlaps = true;
	return res;
}	
