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

/*
 *  Utilities for error reporting
 *
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Authors: Reinhard Hell, Martin Reinecke
 */

#ifndef PLANCK_ERROR_HANDLING_H
#define PLANCK_ERROR_HANDLING_H

#include <string>
#include <iostream>

#if defined (__GNUC__)
#define PLANCK_FUNC_NAME__ __PRETTY_FUNCTION__
#else
#define PLANCK_FUNC_NAME__ 0
#endif

void planck_failure__(const char *file, int line, const char *func,
  const std::string &msg);
void planck_failure__(const char *file, int line, const char *func,
  const char *msg);
void killjob__();

class PlanckError
  {
  private:
    std::string msg;

  public:
    explicit PlanckError(const std::string &message);
    explicit PlanckError(const char *message);

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~PlanckError();
  };

/*! \defgroup errorgroup Error handling */
/*! \{ */

/*! Writes diagnostic output and exits with an error status. */
#define planck_fail(msg) \
do { planck_failure__(__FILE__,__LINE__,PLANCK_FUNC_NAME__,msg); \
throw PlanckError(msg); } while(0)

/*! Throws a PlanckError without diagnostic message. */
#define planck_fail_quietly(msg) \
do { throw PlanckError(msg); } while(0)

/*! Writes diagnostic output and exits with an error status if \a testval
    is \a false. */
#define planck_assert(testval,msg) \
do { if (testval); else planck_fail(msg); } while(0)

/*! Macro for improving error diagnostics. Should be placed immediately
    after the opening brace of \c main(). Must be used in conjunction with
    \c PLANCK_DIAGNOSIS_END. */
#define PLANCK_DIAGNOSIS_BEGIN try {
/*! Macro for improving error diagnostics. Should be placed immediately
    before the closing brace of \c main(). Must be used in conjunction with
    \c PLANCK_DIAGNOSIS_BEGIN. */
#define PLANCK_DIAGNOSIS_END \
} \
catch (PlanckError &) \
  { killjob__(); /* no need for further diagnostics; they were shown already */ } \
catch (std::exception &e) \
  { std::cerr << "std::exception: " << e.what() << std::endl; killjob__(); } \
catch (...) \
  { std::cerr << "Unknown exception" << std::endl; killjob__(); }

/*! \} */

#endif
