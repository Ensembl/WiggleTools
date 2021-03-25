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

#include <stdlib.h> 
#include "wiggletools.h"
#include "hash.h"

const size_t BUCKET_COUNT = 1024;
const size_t BUCKET_INIT_SIZE = 10;

typedef struct hashel_st {
	int key;
	int value;
} HashElement;

typedef struct hashlist_st {
	size_t size;
	size_t count;
	HashElement * list;
} HashList;

struct hash_st {
	HashList * lists;
};

static void hashlist_init(HashList * hashlist) {
	hashlist->size = BUCKET_INIT_SIZE;
	hashlist->list = (HashElement *) calloc(sizeof(HashElement), BUCKET_INIT_SIZE);
}

Hash *hash_construct() {
	Hash * res = (Hash *) calloc(sizeof(Hash), 1);
	res->lists = calloc(sizeof(HashList), BUCKET_COUNT);
	int i;
	for (i =0; i < BUCKET_COUNT; i++)
		hashlist_init(res->lists + i);
	return res;
}

static void hashlist_destroy(HashList * hl) {
	free(hl->list);
}

void hash_destroy(Hash * hash) {
	int i;
	for (i =0; i < BUCKET_COUNT; i++)
		hashlist_destroy(hash->lists + i);
	free(hash->lists);
	free(hash);
}

static bool hashlist_increment(HashList* hl, int key) {
	int i;
	for (i = 0; i < hl->count; i++)
		if (hl->list[i].key == key) {
			hl->list[i].value++;
			return true;
		}

	hl->list[hl->count].key = key;
	hl->list[hl->count].value = 1;
	hl->count++;
	return false;
}

bool hash_increment(Hash * hash, int key) {
	return hashlist_increment(hash->lists + (key % BUCKET_COUNT), key);
}

static int hashlist_remove(HashList * hl, int key) {
	int i;
	int value = 0;
	bool found = false;

	for (i = 0; i < hl->count; i++) {
		if (hl->list[i].key == key) {
			value = hl->list[i].value;
			found = true;
		}
		if (found && i < hl->count - 1) {
			hl->list[i].key = hl->list[i+1].key;
			hl->list[i].value = hl->list[i+1].value;
		}
	}

	if (found) 
		hl->count--;
	return value;
}

int hash_remove(Hash * hash, int key) {
	return hashlist_remove(hash->lists + (key % BUCKET_COUNT), key);
}
	
