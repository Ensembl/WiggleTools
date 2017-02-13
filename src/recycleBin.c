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

/****************************************************************\
*                                                                *
*  Efficient Memory Allocation Routines                          *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2005.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "recycleBin.h"

typedef struct RecycleBin_Node {
	struct RecycleBin_Node *next;
} RecycleBin_Node;

typedef struct chunk_st {
	struct chunk_st *next;
} Chunk;

struct recycleBin_st {
	Chunk *chunk_list;
	RecycleBin_Node *recycle;
	size_t node_size;
	int chunk_pos;
	int nodes_per_chunk;
};

static void initRecycleBin(RecycleBin *recycleBin,
			   size_t node_size, int nodes_per_chunk)
{
	size_t chunckSize, allocSize;

	chunckSize = sizeof(Chunk) + nodes_per_chunk * node_size;
	allocSize = 1;
	/* Get nearest power of 2 */
	while (allocSize < chunckSize)
		allocSize <<= 1;
	nodes_per_chunk = (allocSize - sizeof(Chunk)) / node_size;
	recycleBin->chunk_list = NULL;
	recycleBin->chunk_pos = nodes_per_chunk;
	recycleBin->nodes_per_chunk = nodes_per_chunk;
	recycleBin->node_size = node_size;
	recycleBin->recycle = NULL;
}

RecycleBin *newRecycleBin(size_t node_size, int nodes_per_chunk)
{
	RecycleBin *recycleBin;

	if (node_size < sizeof(RecycleBin_Node)) {
		fprintf(stderr, "Too small elements to create a recycle bin!\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(-1);
	}
	recycleBin = malloc(sizeof(RecycleBin));
	if (!recycleBin) {
		fprintf(stderr, "Out of memory, exiting.\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(-1);
	}
	initRecycleBin (recycleBin, node_size, nodes_per_chunk);

	return recycleBin;
}

static void destroyRecycleBinChunks(RecycleBin * recycleBin)
{
	while (recycleBin->chunk_list != NULL)
	{
		Chunk *chunk;

		chunk = recycleBin->chunk_list;
		recycleBin->chunk_list = recycleBin->chunk_list->next;
		free(chunk);
	}
}

void destroyRecycleBin(RecycleBin * recycleBin)
{
	if (recycleBin == NULL)
		return;

	destroyRecycleBinChunks(recycleBin);
	free(recycleBin);
}

void *allocatePointer(RecycleBin * recycle_bin)
{
	RecycleBin_Node *node;
	Chunk *chunk;

	if (recycle_bin == NULL) {
		fprintf(stderr, "Null recycle bin!\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(-1);
	}

	if (recycle_bin->recycle != NULL) {
		node = recycle_bin->recycle;
		recycle_bin->recycle = node->next;
		return node;
	}

	if (recycle_bin->chunk_pos == recycle_bin->nodes_per_chunk) {
		chunk = malloc(sizeof(Chunk) + recycle_bin->nodes_per_chunk
			       * recycle_bin->node_size);
		if (chunk == NULL) {
			fprintf(stderr, "No more memory for memory chunk!\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(-1);
		}
		chunk->next = recycle_bin->chunk_list;
		recycle_bin->chunk_list = chunk;
		recycle_bin->chunk_pos = 1;
		return (RecycleBin_Node *) ((size_t) (void *) chunk +
					    sizeof(Chunk));
	}

	chunk = recycle_bin->chunk_list;
	return (RecycleBin_Node *) ((size_t) (void *) chunk + sizeof(Chunk)
				    +
				    (recycle_bin->
				     node_size *
				     recycle_bin->chunk_pos++));
}

void deallocatePointer(RecycleBin * recycle_bin, void *data)
{
	RecycleBin_Node *node = data;

	node->next = recycle_bin->recycle;
	recycle_bin->recycle = node;

	return;
}
