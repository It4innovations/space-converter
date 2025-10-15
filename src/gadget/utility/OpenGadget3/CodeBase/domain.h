/*
* @file
* This file is part of the developer version of GADGET3 and contains
* the license conditions for its usage.
*
* @author GADGET-development team, led by Volker Springel and Klaus Dolag.
*
* @section LICENSE
* Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing authors
* (see change logs). All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Received source code may be modified and used as convenient.
*
* 2. Redistributions of source code or in binary form is only possible with
*    explicit agreement of the copyright holders.
*
* 3. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 4. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* 5. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from this
*    software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#ifndef ALLVARS_H
#include "allvars.h"
#endif

void domain_print(void);
void domain_findSplit_work_balanced(int ncpu, int ndomain);
void domain_findSplit_load_balanced(int ncpu, int ndomain);
int domain_sort_loadorigin(const void *a, const void *b);
int domain_sort_segments(const void *a, const void *b);
void domain_assign_load_or_work_balanced(int mode, int multipledomains);
void domain_allocate(void);
int domain_check_memory_bound(int multipledomains);
int domain_compare_key(const void *a, const void *b);
int domain_compare_key(const void *a, const void *b);
int domain_compare_toplist(const void *a, const void *b);
double domain_particle_costfactor(int i);
int domain_countToGo(size_t nlimit);
void domain_Decomposition(int UseAllTimeBins, int SaveKeys);
#ifdef WRITE_KEY_FILES
int domain_decompose(int SaveKeys);
int domain_determineTopTree(int SaveKeys, int multipledomains);
void domain_findExtent(int SaveKeys);
#else
int domain_decompose(void);
int domain_determineTopTree(void);
void domain_findExtent(void);
#endif
void domain_exchange(void);
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void domain_findSplit(int cpustart, int ncpu, int first, int last);
void domain_findSplit_balanced(int cpustart, int ncpu, int first, int last);
void domain_free(void);
void domain_shiftSplit(void);
void domain_sumCost(void);
int domain_topsplit(int node, peanokey startkey);
int domain_topsplit_local(int node, peanokey startkey, int mode);
int domain_topsplit_special(void);
int domain_compare_key(const void *a, const void *b);
#ifdef AR_PARALLEL_LOCAL_REFINE
int domain_check_for_local_refine(int i, double countlimit, double costlimit, int layer);
#else
int domain_check_for_local_refine(int i, double countlimit, double costlimit);
#endif
void domain_free_trick(void);
void domain_allocate_trick(void);
int domain_recursively_combine_topTree(int start, int ncpu);
void domain_walktoptree(int no);
void mysort_domain(void *b, size_t n, size_t s);

#ifdef RECOMPOSE_DOMAIN
void domain_recomposition();
#endif

#endif //#ifndef __DOMAIN_H__
