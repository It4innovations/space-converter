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
#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdlib.h>
#include <stdio.h>

#ifndef GADGET3_IO_LIB
#include "../gadgetconfig.h"
#endif

#ifdef PARALLEL
int ThisTask;
int NTask;
#endif

#define safe_fgets( str, num, stream ) util_fgets( str, num, stream, __FILE__, __LINE__ )
#define safe_fread( ptr, size, count, stream ) util_fread( ptr, size, count, stream, __FILE__, __LINE__ )

void myprintf(const char *format, ...);
char *util_fgets(char *str, int num, FILE * stream, char *file, int line);
size_t util_fread(void *ptr, size_t size, size_t count, FILE * stream, char *file, int line);

__inline static int imin(int a, int b)
{
  return a < b ? a : b;
}

__inline static int imax(int a, int b)
{
  return a > b ? a : b;
}

__inline static double dmin(double a, double b)
{
  return a < b ? a : b;
}

__inline static double dmax(double a, double b)
{
  return a > b ? a : b;
}

#define VECT_NORM(N, V) ({double VN; int vcc; for(vcc = 0, VN = 0; vcc < (N); vcc++) VN+=((V)[vcc]*(V)[vcc]); VN;})



#endif /* UTILITIES_H */
