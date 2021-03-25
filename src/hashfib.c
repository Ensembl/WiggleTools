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
#include "fib.h"
#include "hashfib.h"

struct hashfib_st {
	Hash *hash;
	FibHeap *fib;
};

HashFib *hashfib_construct() {
	HashFib * res = (HashFib*) calloc(sizeof(HashFib), 1);
	res->hash = hash_construct(); 
	res->fib = fh_makeheap();
	return res;
}

void hashfib_insert(HashFib * hf, int key) {
	if (!hash_increment(hf->hash, key))
		fh_insert(hf->fib, key, 0);
}

bool hashfib_empty(HashFib * hf) {
	return fh_empty(hf->fib);
}

int hashfib_min(HashFib * hf) {
	return fh_min(hf->fib);
}

int hashfib_remove_min(HashFib * hf) {
	int key = fh_min(hf->fib);
	fh_extractmin(hf->fib);
	return hash_remove(hf->hash, key); 
}

void hashfib_destroy(HashFib * hf) {
	hash_destroy(hf->hash);
	fh_deleteheap(hf->fib);
	free(hf);
}
