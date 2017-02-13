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
 * Copyright 1997-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
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
 *	$Id: fib.c,v 1.10 2007/10/19 13:09:26 zerbino Exp $
 *
 */
#include <limits.h>
#include <stdlib.h>

#include "fib.h"
#include "recycleBin.h"

#include "fibpriv.h"

#define BLOCKSIZE 10000

static FibHeapNode *allocateFibHeapEl(FibHeap * heap)
{
	return allocatePointer(heap->nodeMemory);
}

static void deallocateFibHeapEl(FibHeapNode * a, FibHeap * heap)
{
	deallocatePointer(heap->nodeMemory, a);
}

#define INT_BITS        (sizeof(int) * 8)

static inline int ceillog2(int a)
{
	int oa;
	int i;
	int b;
	int cons;

	oa = a;
	b = INT_BITS / 2;
	i = 0;
	while (b) {
		i = (i << 1);
		cons = ((int) 1) << b;
		if (a >= cons) {
			a /= cons;
			i = i | 1;
		} else
			a &= cons - 1;
		b /= 2;
	}
	if ((((int) 1 << i)) == oa)
		return i;
	else
		return i + 1;
}

/*
 * Public Heap Functions
 */
FibHeap *fh_makeheap()
{
	FibHeap *new = malloc(sizeof(FibHeap));

	if (new) {
		new->nodeMemory = newRecycleBin(sizeof(FibHeapNode), BLOCKSIZE);
		new->fh_neginf = NULL;
		new->fh_n = 0;
		new->fh_Dl = -1;
		new->fh_cons = NULL;
		new->fh_min = NULL;
		new->fh_root = NULL;
	}

	return new;
}

void fh_deleteheap(FibHeap * h)
{
	destroyRecycleBin(h->nodeMemory);
	h->fh_neginf = NULL;
	if (h->fh_cons != NULL)
		free(h->fh_cons);
	h->fh_cons = NULL;
	free(h);
}

/*
 * Public Key Heap Functions
 */
FibHeapNode *fh_insert(FibHeap * h, int key, int value)
{
	FibHeapNode *x;

	if ((x = fhe_newelem(h)) == NULL)
		return NULL;

	/* just insert on root list, and make sure it's not the new min */
	x->fhe_key = key;
	x->fhe_value = value;

	fh_insertel(h, x);

	return x;
}

static void fh_insertel(FibHeap * h, FibHeapNode * x)
{
	fh_insertrootlist(h, x);

	if (h->fh_min == NULL || x->fhe_key < h->fh_min->fhe_key)
		h->fh_min = x;

	h->fh_n++;
}

static void fh_insertrootlist(FibHeap * h, FibHeapNode * x)
{
	if (h->fh_root == NULL) {
		h->fh_root = x;
		x->fhe_left = x;
		x->fhe_right = x;
	} else {
		fhe_insertafter(h->fh_root, x);
	}
}

static void fhe_insertafter(FibHeapNode * a, FibHeapNode * b)
{
	if (a == a->fhe_right) {
		a->fhe_right = b;
		a->fhe_left = b;
		b->fhe_right = a;
		b->fhe_left = a;
	} else {
		b->fhe_right = a->fhe_right;
		a->fhe_right->fhe_left = b;
		a->fhe_right = b;
		b->fhe_left = a;
	}
}

int fh_min(FibHeap * h)
{
	if (h->fh_min == NULL)
		return (int) INT_MIN;
	return h->fh_min->fhe_key;
}

int fh_notempty(FibHeap * h)
{
	return (int) (h->fh_min != NULL);
}

int fh_empty(FibHeap * h)
{
	return (int) (h->fh_min == NULL);
}

int fh_extractmin(FibHeap * h)
{
	if (h->fh_min != NULL) {
		FibHeapNode * min = fh_extractminel(h);
		int res = min->fhe_value;
		deallocateFibHeapEl(min, h);
		return res;
	}

	return -1;
}

static FibHeapNode *fh_extractminel(FibHeap * h)
{
	FibHeapNode *ret;
	FibHeapNode *x, *y, *orig;

	ret = h->fh_min;

	orig = NULL;
	/* put all the children on the root list */
	/* for true consistancy, we should use fhe_remove */
	for (x = ret->fhe_child; x != orig && x != NULL;) {
		if (orig == NULL)
			orig = x;
		y = x->fhe_right;
		x->fhe_p = NULL;
		fh_insertrootlist(h, x);
		x = y;
	}
	/* remove minimum from root list */
	fh_removerootlist(h, ret);
	h->fh_n--;

	/* if we aren't empty, consolidate the heap */
	if (h->fh_n == 0)
		h->fh_min = NULL;
	else {
		h->fh_min = ret->fhe_right;
		fh_consolidate(h);
	}

	return ret;
}

static void fh_removerootlist(FibHeap * h, FibHeapNode * x)
{
	if (x->fhe_left == x)
		h->fh_root = NULL;
	else
		h->fh_root = fhe_remove(x);
}

static void fh_consolidate(FibHeap * h)
{
	FibHeapNode **a;
	FibHeapNode *w;
	FibHeapNode *y;
	FibHeapNode *x;
	int i;
	int d;
	int D;

	fh_checkcons(h);

	/* assign a the value of h->fh_cons so I don't have to rewrite code */
	D = h->fh_Dl + 1;
	a = h->fh_cons;

	for (i = 0; i < D; i++)
		a[i] = NULL;

	while ((w = h->fh_root) != NULL) {
		x = w;
		fh_removerootlist(h, w);
		d = x->fhe_degree;
		/* XXX - assert that d < D */
		while (a[d] != NULL) {
			y = a[d];
			if (fh_compare(h, x, y) > 0) {
				FibHeapNode * temp = x;
				x = y;
				y = temp;
			}
			fh_heaplink(h, y, x);
			a[d] = NULL;
			d++;
		}
		a[d] = x;
	}
	h->fh_min = NULL;
	for (i = 0; i < D; i++)
		if (a[i] != NULL) {
			fh_insertrootlist(h, a[i]);
			if (h->fh_min == NULL
			    || fh_compare(h, a[i], h->fh_min) < 0)
				h->fh_min = a[i];
		}
}

static void fh_checkcons(FibHeap * h)
{
	int oDl;

	/* make sure we have enough memory allocated to "reorganize" */
	if (h->fh_Dl == -1 || h->fh_n > (1 << h->fh_Dl)) {
		oDl = h->fh_Dl;
		if ((h->fh_Dl = ceillog2(h->fh_n) + 1) < 8)
			h->fh_Dl = 8;
		if (oDl != h->fh_Dl)
			h->fh_cons =
			    (FibHeapNode **) realloc(h->fh_cons,
						     sizeof *h->
						     fh_cons *
						     (h->fh_Dl + 1));
		if (h->fh_cons == NULL)
			abort();
	}
}

static int fh_compare(FibHeap * h, FibHeapNode * a, FibHeapNode * b)
{
	if (a->fhe_key < b->fhe_key)
		return -1;
	if (a->fhe_key == b->fhe_key)
		return 0;
	return 1;
}

static void fh_heaplink(FibHeap * h, FibHeapNode * y, FibHeapNode * x)
{
	/* make y a child of x */
	if (x->fhe_child == NULL)
		x->fhe_child = y;
	else
		fhe_insertafter(x->fhe_child->fhe_left, y);
	y->fhe_p = x;
	x->fhe_degree++;
	y->fhe_mark = 0;
}

/*
 * begining of handling elements of fibheap
 */
static FibHeapNode *fhe_newelem(FibHeap * h)
{
	FibHeapNode *e;

	if ((e = allocateFibHeapEl(h)) == NULL)
		return NULL;

	e->fhe_degree = 0;
	e->fhe_mark = 0;
	e->fhe_p = NULL;
	e->fhe_child = NULL;
	e->fhe_left = e;
	e->fhe_right = e;

	return e;
}

static FibHeapNode *fhe_remove(FibHeapNode * x)
{
	FibHeapNode *ret;

	if (x == x->fhe_left)
		ret = NULL;
	else
		ret = x->fhe_left;

	/* fix the parent pointer */
	if (x->fhe_p != NULL && x->fhe_p->fhe_child == x)
		x->fhe_p->fhe_child = ret;

	x->fhe_right->fhe_left = x->fhe_left;
	x->fhe_left->fhe_right = x->fhe_right;

	/* clear out hanging pointers */
	x->fhe_p = NULL;
	x->fhe_left = x;
	x->fhe_right = x;

	return ret;
}
