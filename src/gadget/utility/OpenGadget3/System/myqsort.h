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
/* must define macros QSORT, KEY_TYPE, STRUCT_TYPE, KEY_COPY, GET_KEYVAL */
static void QSORT(KEY_TYPE *a,int n)
{
#pragma alloca
  KEY_TYPE *b = (KEY_TYPE *)alloca(sizeof(KEY_TYPE)*n);

  if (n > 1) {
     int i,j1,j2;
     double avg = 0.0;
     KEY_BASE_TYPE kavg;

     /* compute average key value and use as pivot point */
     
     for(i=0;i < n;i++)
        avg += KEY_GETVAL(&a[i]);

     kavg = (KEY_BASE_TYPE) (avg / (double)n);

     j1 = 0;
     j2 = n;

     for (i = 0;i < n;i++)
       if (KEY_GETVAL(&a[i]) <= kavg) {
          KEY_COPY(&a[i],&b[j1]);
          j1++;
       }
       else {
          --j2;
          KEY_COPY(&a[i],&b[j2]);
       }

     if (j1 > 0 && j2 < n) {
     
        if (j1 > 0) QSORT(b,j1);
        if (j2 < n) QSORT(&b[j1],n-j1);

        memcpy(a,b,sizeof(KEY_TYPE)*n);
     }
  }
}
