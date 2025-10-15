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
#ifndef _BIGFLOAT_H
#define _BIGFLOAT_H

#if (GDE_BIGFLOAT == 0)
#include <xpre.h>

class BigFloat 
{
public:
  struct xpr val;

  /* self + */
  BigFloat & operator+=(const BigFloat & h)
  {
    val = xadd(val, h.val, 0);
    return *this;
  };

  BigFloat & operator+=(double d)
  {
    val = xadd(val, dbltox(d), 0);
    return *this;
  };

  BigFloat & operator+=(float f)
  {
    val = xadd(val, flttox(f), 0);
    return *this;
  };

  /* self - */
  BigFloat & operator-=(const BigFloat & h)
  {
    val = xadd(val, h.val, 1);
    return *this;
  };

  BigFloat & operator-=(double d)
  {
    val = xadd(val, dbltox(d), 1);
    return *this;
  };

  BigFloat & operator-=(float f)
  {
    val = xadd(val, flttox(f), 1);
    return *this;
  };

  /* self * */
  BigFloat & operator*=(const BigFloat & h)
  {
    val = xmul(val, h.val);
    return *this;
  };

  BigFloat & operator*=(double d)
  {
    val = xmul(val, dbltox(d));
    return *this;
  };

  BigFloat & operator*=(float f)
  {
    val = xmul(val, flttox(f));
    return *this;
  };

  /* self / */
  BigFloat & operator/=(const BigFloat & h)
  {
    val = xdiv(val, h.val);
    return *this;
  };

  BigFloat & operator/=(double d)
  {
    val = xdiv(val, dbltox(d));
    return *this;
  };

  BigFloat & operator/=(float f)
  {
    val = xdiv(val, flttox(f));
    return *this;
  };

  /* + */
  BigFloat operator+(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = xadd(val, h.val, 0);
    return htmp;
  }

  BigFloat operator+(double d)
  {
    BigFloat htmp;
    htmp.val = xadd(val, dbltox(d), 0);
    return htmp;
  }

  BigFloat operator+(float f)
  {
    BigFloat htmp;
    htmp.val = xadd(val, flttox(f), 0);
    return htmp;
  }

  /* - */
  BigFloat operator-(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = xadd(val, h.val, 1);
    return htmp;      
  }

  BigFloat operator-(double d)
  {
    BigFloat htmp;
    htmp.val = xadd(val, dbltox(d), 1);
    return htmp;
  }

  BigFloat operator-(float f)
  {
    BigFloat htmp;
    htmp.val = xadd(val, flttox(f), 1);
    return htmp;
  }

  /* * */
  BigFloat operator*(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = xmul(val, h.val);
    return htmp;      
  }

  BigFloat operator*(double d)
  {
    BigFloat htmp;
    htmp.val = xmul(val, dbltox(d));
    return htmp;
  }

  BigFloat operator*(float f)
  {
    BigFloat htmp;
    htmp.val = xmul(val, flttox(f));
    return htmp;
  }

  /* / */
  BigFloat operator/(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = xdiv(val, h.val);
    return htmp; 
  }

  BigFloat operator/(double d)
  {
    BigFloat htmp;
    htmp.val = xdiv(val, dbltox(d));
    return htmp;
  }

  BigFloat operator/(float f)
  {
    BigFloat htmp;
    htmp.val = xdiv(val, flttox(f));
    return htmp;
  }

  BigFloat abs(void)
  {
   BigFloat htmp;
   htmp.val = xabs(val);
   return htmp;
  }

  BigFloat log(void)
  {
   BigFloat htmp;
   htmp.val = xlog(val);
   return htmp;
  }
 
  BigFloat sqrt()
  { 
   BigFloat htmp;
   htmp.val = xsqrt(val);
   return htmp;
  }

  /* comparison */
  bool operator==(const BigFloat & h)
   {
    return (xprcmp(&val, &h.val)==0)?true:false;
   }

  bool operator<(const BigFloat & h)
   {
    return (xprcmp(&val, &h.val)==-1)?true:false;
   }

  bool operator>(const BigFloat & h)
   {
    return (xprcmp(&val, &h.val)==+1)?true:false;
   }

  bool operator<=(const BigFloat & h)
   {
    return (xgt(val, h.val)==0)?true:false;
   }

  bool operator>=(const BigFloat & h)
   {
    return (xlt(val, h.val)==0)?true:false;
   }

  bool operator!=(const BigFloat & h)
   {
    return (xeq(val, h.val)==0)?true:false;
   }



  /* assign */
  BigFloat & operator=(const BigFloat & h)
    {
      val = h.val;
      return *this;
    };

  BigFloat & operator=(double d)
    {
      val = dbltox(d);
      return *this;
    };

  BigFloat & operator=(float f)
    {
      val = flttox(f);
      return *this;
    };




  /* cast */
  operator  double () const 
   {
    return xtodbl(val);
   };

  operator  float () const
   {
    return (float)xtodbl(val);
   };

 

};
#endif

#if (GDE_BIGFLOAT == 1)
#include "ttmath/ttmath.h"

//ttmath:Big<words for exponent, words for mantissa>
//typedef ttmath::Big<1,15> MyBig;
typedef ttmath::Big<1,10> MyBig;

class BigFloat 
{
public:
  MyBig val;


  /* self + */
  BigFloat & operator+=(const BigFloat & h)
  {
    val.Add(h.val);
    return *this;
  };

  BigFloat & operator+=(double d)
  {
    val.Add(d);
    return *this;
  };

  BigFloat & operator+=(float f)
  {
    val.Add(f);
    return *this;
  };

  /* self - */
  BigFloat & operator-=(const BigFloat & h)
  {
    val.Sub(h.val);
    return *this;
  };

  BigFloat & operator-=(double d)
  {
    val.Sub(d);
    return *this;
  };

  BigFloat & operator-=(float f)
  {
    val.Sub(f);
    return *this;
  };

  /* self * */
  BigFloat & operator*=(const BigFloat & h)
  {
    val.Mul(h.val);
    return *this;
  };

  BigFloat & operator*=(double d)
  {
    val.Mul(d);
    return *this;
  };

  BigFloat & operator*=(float f)
  {
    val.Mul(f);
    return *this;
  };

  /* self / */
  BigFloat & operator/=(const BigFloat & h)
  {
    val.Div(h.val);
    return *this;
  };
    
  BigFloat & operator/=(double d)
  {
    val.Div(d);
    return *this;
  };

  BigFloat & operator/=(float f)
  {
    val.Div(f);
    return *this;
  };

  /* + */
  BigFloat operator+(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Add(h.val);
    return htmp;
  }

  BigFloat operator+(double d)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Add(d);
    return htmp;
  }

  BigFloat operator+(float f)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Add(f);
    return htmp;
  }

  /* - */
  BigFloat operator-(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Sub(h.val);
    return htmp;
  }

  BigFloat operator-(double d)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Sub(d);
    return htmp;
  }

  BigFloat operator-(float f)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Sub(f);
    return htmp;
  }

  /* * */
  BigFloat operator*(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Mul(h.val);
    return htmp;
  }

  BigFloat operator*(double d)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Mul(d);
    return htmp;
  }

  BigFloat operator*(float f)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Mul(f);
    return htmp;
  }

  /* / */
  BigFloat operator/(const BigFloat & h)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Div(h.val);
    return htmp;
  }

  BigFloat operator/(double d)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Div(d);
    return htmp;
  }

  BigFloat operator/(float f)
  {
    BigFloat htmp;
    htmp.val = val;
    htmp.val.Div(f);
    return htmp;
  }

  BigFloat abs(void)
  {
   BigFloat htmp;
   htmp.val = Abs(val);
   return htmp;
  }

  BigFloat log(void)
  {
   BigFloat htmp;
   htmp.val = Ln(val);
   return htmp;

  }

  BigFloat sqrt()
  {
   BigFloat htmp;
   htmp.val = Sqrt(val);
   return htmp;
  }

  /* comparison */
  bool operator==(const BigFloat & h)
   {
    return (val==h.val);
   }

  bool operator<(const BigFloat & h)
   {
    return (val<h.val); 
   }

  bool operator>(const BigFloat & h)
   {
    return (val>h.val);
   }

  bool operator<=(const BigFloat & h)
   {
    return (val<=h.val);
   }

  bool operator>=(const BigFloat & h)
   {
    return (val>=h.val);
   }

  bool operator!=(const BigFloat & h)
   {
    return (val!=h.val);
   }


  /* assign */
  BigFloat & operator=(const BigFloat & h)
    {
      val = h.val;
      return *this;
    };

  BigFloat & operator=(double d)
    {
      val = d;
      return *this;
    };

  BigFloat & operator=(float f)
    {
      val = f;
      return *this;
    };




  /* cast */
  operator  double () const
   {
    return val.ToDouble();
   };

  operator  float () const
   {
    return val.ToFloat();
   };
};


#endif
#endif 
