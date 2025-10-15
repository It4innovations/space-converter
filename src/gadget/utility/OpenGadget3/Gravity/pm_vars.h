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

#ifndef GRIDBOOST
#define GRIDBOOST 2
#endif
#define  GRID  (GRIDBOOST*PMGRID)
#define  GRID2 (2*(GRID/2 + 1))
#define  PMGRID2 (2*(PMGRID/2 + 1))


#if (GRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif


#define d_fftw_real fftw_real

extern rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

extern int slab_to_task[GRID];
extern int *slabs_per_task;
extern int *first_slab_of_task;

extern int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

extern int fftsize, maxfftsize;

extern fftw_real *kernel[2], *rhogrid, *forcegrid, *workspace;
extern fftw_complex *fft_of_kernel[2], *fft_of_rhogrid;
extern d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;

extern int *part_sortindex;

extern struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} *part;
