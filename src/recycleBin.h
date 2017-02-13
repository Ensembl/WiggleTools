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

#ifndef INCLUDED_RECYCLEBIN_H
#define INCLUDED_RECYCLEBIN_H

typedef struct recycleBin_st RecycleBin;

// Constructor, Destructor
RecycleBin *newRecycleBin(size_t node_size, int nodes_per_chunk);
void destroyRecycleBin(RecycleBin * recycle_bin);

// Use
void *allocatePointer(RecycleBin * recycle_bin);
void deallocatePointer(RecycleBin * recycle_bin, void *data);

#endif				/* INCLUDED_RECYCLEBIN_H */
