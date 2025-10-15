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
/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file safe_cast.h
 *  Numerical cast operator with additional checks that the value is preserved.
 *
 *  Copyright (C) 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_SAFE_CAST_H
#define PLANCK_SAFE_CAST_H

#include <limits>
#include "error_handling.h"

template<typename T1, typename T2, bool s1, bool s2> struct safe_cast_helper__
  {};

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,true,true>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert(T2(res)==arg, "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,false,false>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert(T2(res)==arg, "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,true,false>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert((res>=0) && (T2(res)==arg),
      "safe_cast: value changed during cast");
    return res;
    }
  };

template<typename T1, typename T2> struct safe_cast_helper__ <T1,T2,false,true>
  {
  static T1 cast (const T2 &arg)
    {
    T1 res = T1(arg);
    planck_assert((arg>=0) && (T2(res)==arg),
      "safe_cast: value changed during cast");
    return res;
    }
  };

/*! Tries to cast \a arg from its type to a variable of type \c T1.
    If this conversion leads to a change of the actual value (e.g. due to
    overflow or truncation), an exception is thrown. */
template<typename T1, typename T2> inline T1 safe_cast(const T2 &arg)
  {
  return safe_cast_helper__<T1,T2,std::numeric_limits<T1>::is_signed,
    std::numeric_limits<T2>::is_signed>::cast(arg);
  }

#endif
