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
#ifndef _VECTOR_H_
#define _VECTOR_H_

#if defined(VECTOR_SSE)
#include "vector_sse2.h"
#elif defined(VECTOR_AVX)
#include "vector_avx.h"
#elif defined(VECTOR_VSX)
#include "vector_vsx.h"
#elif defined(VECTOR_QPX)
#include "vector_qpx.h"
#else

// generic scalar code

#ifdef DOUBLEPRECISION

typedef ALIGN(32) union {
  unsigned long long i[4];
  double d[4];
} t_vector;

#else

typedef ALIGN(16) union {
  unsigned int i[4];
  float d[4];
} t_vector;

#endif

// read the first 3 components of a vector, 4th component is undefined
static inline void LOAD_VECTOR3(MyLongDouble * src,t_vector * v)
{
  v->d[0] = src[0];
  v->d[1] = src[1];
  v->d[2] = src[2];
}

// read all 4 components of a vector
static inline void LOAD_VECTOR4(MyLongDouble * src,t_vector * v)
{
  v->d[0] = src[0];
  v->d[1] = src[1];
  v->d[2] = src[2];
  v->d[3] = src[3];
}

// stores only the first 3 components of a vector (avoid, because it's slow)
static inline void STORE_VECTOR3(MyLongDouble * dst,t_vector * v)
{
  dst[0] = v->d[0];
  dst[1] = v->d[1];
  dst[2] = v->d[2];
}

// stores all 4 components of a vector (fast)
static inline void STORE_VECTOR4(MyLongDouble * dst,t_vector * v)
{
  dst[0] = v->d[0];
  dst[1] = v->d[1];
  dst[2] = v->d[2];
  dst[3] = v->d[3];
}

// initializes first 3 components of a vector (4th component becomes undefined)
static inline void SET_VECTOR3(MyLongDouble a,t_vector * v)
{
  v->d[0] = a;
  v->d[1] = a;
  v->d[2] = a;
}

// initializes all 4 components of a vector with the same scalar value
static inline void SET_VECTOR4(MyLongDouble a,t_vector * v)
{
  v->d[0] = a;
  v->d[1] = a;
  v->d[2] = a;
  v->d[3] = a;
}

// initializes first 3 components of a vector (4th component becomes undefined)
static inline void INIT_VECTOR3(MyLongDouble a0,MyLongDouble a1,MyLongDouble a2,t_vector * v)
{
  v->d[0] = a0;
  v->d[1] = a1;
  v->d[2] = a2;
}

// adds the first 3 components of 2 vectors (4th component undefined)
static inline void ADD_VECTOR3(t_vector * a,t_vector * b,t_vector * result)
{
  result->d[0] = a->d[0]+b->d[0];
  result->d[1] = a->d[1]+b->d[1];
  result->d[2] = a->d[2]+b->d[2];
}

// multiplies the first 3 components of 2 vectors (4th component undefined)
static inline void MUL_VECTOR3(t_vector * a,t_vector * b,t_vector * result)
{
  result->d[0] = a->d[0] * b->d[0];
  result->d[1] = a->d[1] * b->d[1];
  result->d[2] = a->d[2] * b->d[2];
}

// multiplies the first 3 components of 2 vectors (4th component undefined)
static inline void SCALE_VECTOR3(MyLongDouble a,t_vector * b,t_vector * v)
{
  v->d[0] = a * b->d[0];
  v->d[1] = a * b->d[1];
  v->d[2] = a * b->d[2];
}

// returns a value >0, if a(i) < b(i) for any i=1..3
static inline int ANY_COMP_LT_VECTOR3(t_vector * a,t_vector * b,int mask)
{
  return (a->d[0] < b->d[0]) | (a->d[1] < b->d[1]) | (a->d[2] < b->d[2]);
}

// returns 0, if a(i) < b(i) for any i=1..3
static inline int ALL_COMP_LT_VECTOR3(t_vector * a,t_vector * b,int mask)
{
  return (a->d[0] < b->d[0]) & (a->d[1] < b->d[1]) & (a->d[2] < b->d[2]);
}

// returns the L2 norm of the first 3 components of a vector
static inline MyDouble L2NORM_VECTOR3(t_vector * v)
{
  return v->d[0] * v->d[0] + v->d[1] * v->d[1] + v->d[2] * v->d[2];
}

#endif

#ifdef DO_NOT_BRACH_IF

// on the Power platform it is more efficient to compute all conditions first
// and then perform a single branch
// On current Intel processors this is equally fast

inline static int ngb_check_node(struct NODE * cur,t_vector * v2,t_vector * box,t_vector * hbox,MyDouble hsml)
{
  MyDouble dist  = hsml + 0.5*cur->len;
  MyDouble dist2 = dist + FACT1*cur->len;
  MyDouble d2,dx,dy,dz;
  int node = cur->u.d.sibling; // in case the node can be discarded

  dx = NGB_PERIODIC_LONG(cur->center[0] - v2->d[0],box->d[0],hbox->d[0]);
  dy = NGB_PERIODIC_LONG(cur->center[1] - v2->d[1],box->d[1],hbox->d[1]);
  dz = NGB_PERIODIC_LONG(cur->center[2] - v2->d[2],box->d[2],hbox->d[2]);
  d2 = dx * dx + dy * dy + dz * dz;
  // now test against the minimal sphere enclosing everything
  if ((dx <= dist) & (dy <= dist) & (dz <= dist) & (d2 <= dist2 * dist2)) node = cur->u.d.nextnode;     // ok, we need to open the node
  return node;
}

#else

// On older Intel and AMD processors it seems better to avoid the computations and branch early
// TODO: This should be validated

inline static int ngb_check_node(struct NODE * cur,t_vector * v2,t_vector * box,t_vector * hbox,MyDouble hsml)
{
  MyLongDouble dist  = hsml + 0.5*cur->len;
  MyLongDouble d2,dx,dy,dz;
  int node = cur->u.d.sibling; // in case the node can be discarded

  dx = NGB_PERIODIC_LONG(cur->center[0] - v2->d[0],box->d[0],hbox->d[0]);
  if (dx <= dist) {
    dy = NGB_PERIODIC_LONG(cur->center[1] - v2->d[1],box->d[1],hbox->d[1]);
    if (dy <= dist) {
      dz = NGB_PERIODIC_LONG(cur->center[2] - v2->d[2],box->d[2],hbox->d[2]);
      if (dz <= dist) {
        // now test against the minimal sphere enclosing everything
        dist += FACT1 * cur->len;
        d2 = dx * dx + dy * dy + dz * dz;
        if (d2 <= dist * dist) node = cur->u.d.nextnode;     // ok, we need to open the node
      }
    }
  }
  return node;
}

#endif

#endif
