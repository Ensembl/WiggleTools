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

/*-
 * Copyright 1997, 1999-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without

 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fibpriv.h,v 1.10 2007/10/09 09:56:46 zerbino Exp $
 *
 */

#ifndef _FIBPRIV_H_
#define _FIBPRIV_H_

#ifndef bool
#define bool char
#define true 1
#define false 0
#endif

/*
 * specific node operations
 */
struct fibheap_el {
	FibHeapNode *fhe_p;
	FibHeapNode *fhe_child;
	FibHeapNode *fhe_left;
	FibHeapNode *fhe_right;
	int fhe_key;
	int fhe_value;
	int fhe_degree;
	bool fhe_mark;
};

static FibHeapNode *fhe_newelem(struct fibheap *);
static void fhe_insertafter(FibHeapNode * a, FibHeapNode * b);
static FibHeapNode *fhe_remove(FibHeapNode * a);

/*
 * global heap operations
 */
struct fibheap {
	RecycleBin *nodeMemory;
	int fh_n;
	int fh_Dl;
	FibHeapNode **fh_cons;
	FibHeapNode *fh_min;
	FibHeapNode *fh_root;
	void *fh_neginf;
	bool fh_keys;
};

static void fh_insertrootlist(FibHeap * h, FibHeapNode * x);
static void fh_removerootlist(FibHeap *, FibHeapNode *);
static void fh_consolidate(FibHeap *);
static void fh_heaplink(FibHeap * h, FibHeapNode * y, FibHeapNode * x);
static FibHeapNode *fh_extractminel(FibHeap *);
static void fh_checkcons(FibHeap * h);
static int fh_compare(FibHeap * h, FibHeapNode * a, FibHeapNode * b);
static void fh_insertel(FibHeap * h, FibHeapNode * x);

/*
 * general functions
 */
static inline int ceillog2(int a);

#endif				/* _FIBPRIV_H_ */
