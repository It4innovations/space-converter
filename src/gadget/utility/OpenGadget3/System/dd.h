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

/* I (Volker Springel) have written this class `dd' based in part on 
 * code contained in the 
 *  ------------------------------------------------------------------------
 *  | QUAD-DOUBLE/DOUBLE-DOUBLE COMPUTATION PACKAGE                        |
 *  |                                                                      |
 *  | Yozo Hida        U.C. Berkeley               yozo@cs.berkeley.edu    |
 *  | Xiaoye S. Li     Lawrence Berkeley Natl Lab  xiaoye@nersc.gov        |
 *  | David H. Bailey  Lawrence Berkeley Natl Lab  dhbailey@lbl.gov        |
 *  |                                                                      |
 *  | Revised  2005-03-12  Copyright (c) 2005                              |
 *  ------------------------------------------------------------------------
 *  (available at http://crd.lbl.gov/~dhbailey/mpdist/)
 */

#ifndef _QD_DD_H
#define _QD_DD_H

class dd
{
public:
  double hi, lo;

  double quick_two_sum(double a, double b, double &err)
  {
    double s = a + b;

    err = b - (s - a);
    return s;
  }

  double two_sum(double a, double b, double &err)
  {
    double s = a + b;
    double bb = s - a;

    err = (a - (s - bb)) + (b - bb);
    return s;
  }


  /* Self-Addition with a double */
  dd & operator+=(double a)
  {
    double s1, s2;

    s1 = two_sum(hi, a, s2);
    s2 += lo;
    hi = quick_two_sum(s1, s2, lo);
    return *this;
  };

  /* Self-Addition with a doubledouble */
  dd & operator+=(const dd & a)
  {
    double s1, s2, t1, t2;

    s1 = two_sum(hi, a.hi, s2);
    t1 = two_sum(lo, a.lo, t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, s2);
    s2 += t2;
    hi = quick_two_sum(s1, s2, lo);
    return *this;
  };

  /* Assignment */
  dd & operator=(double a)
    {
      hi = a;
      lo = 0.0;
      return *this;
    };

  /* Cast */
  operator  double () const 
   {
    return hi;
   };

};

#endif /* _QD_DD_H */
