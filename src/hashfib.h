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

#ifndef _HASHFIB_H_
#define _HASHFIB_H_

typedef struct hashfib_st HashFib;

HashFib *hashfib_construct();
bool hashfib_empty(HashFib*);
void hashfib_insert(HashFib *, int);
int hashfib_min(HashFib *);
int hashfib_remove_min(HashFib *);
void hashfib_destroy(HashFib *);

#endif				/* _FIB_H_ */
