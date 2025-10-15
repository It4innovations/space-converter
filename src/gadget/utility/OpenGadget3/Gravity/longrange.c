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
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"


/*! \file longrange.c
 *  \brief driver routines for computation of long-range gravitational PM force
 */

#ifdef PMGRID

/*! Driver routine to call initializiation of periodic or/and non-periodic FFT
 *  routines.
 */
void long_range_init(void)
{
  pm_init();
}


void long_range_init_regionsize(void)
{

#ifdef PERIODIC
  if(RestartFlag != 1)
    pm_init_regionsize();
#ifdef PLACEHIGHRESREGION
  pm_setup_nonperiodic_kernel();
#endif
#else
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif

}


/*! This function computes the long-range PM force for all particles.
 */
void long_range_force(void)
{
  int i;

#ifndef PERIODIC
  int j;
  double fac;
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif
#ifdef DISTORTIONTENSORPS
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif
    }

#ifdef NOGRAVITY
  return;
#endif


#ifdef PERIODIC
  i = pmforce(2);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
  pmtidaltensor_periodic_fourier(0);
  pmtidaltensor_periodic_fourier(1);
  pmtidaltensor_periodic_fourier(2);
  pmtidaltensor_periodic_fourier(3);
  pmtidaltensor_periodic_fourier(4);
  pmtidaltensor_periodic_fourier(5);
*/
/* FINITE DIFFERENCES based */
  pmtidaltensor_periodic_diff();
  /* check tidal tensor of given particle ID */
//  check_tidaltensor_periodic(10000);
#endif
#ifdef PLACEHIGHRESREGION
  i = pmforce(1);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
  pmtidaltensor_nonperiodic_fourier(1, 0);
  pmtidaltensor_nonperiodic_fourier(1, 1);
  pmtidaltensor_nonperiodic_fourier(1, 2);
  pmtidaltensor_nonperiodic_fourier(1, 3);
  pmtidaltensor_nonperiodic_fourier(1, 4);
  pmtidaltensor_nonperiodic_fourier(1, 5);
*/
/* FINITE DIFFERENCES based */
  pmtidaltensor_nonperiodic_diff(1);
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      /* try again */
      i = pmforce(1);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
      pmtidaltensor_nonperiodic_fourier(1, 0);
      pmtidaltensor_nonperiodic_fourier(1, 1);
      pmtidaltensor_nonperiodic_fourier(1, 2);
      pmtidaltensor_nonperiodic_fourier(1, 3);
      pmtidaltensor_nonperiodic_fourier(1, 4);
      pmtidaltensor_nonperiodic_fourier(1, 5);
*/
/* FINITE DIFFERENCES based */
      pmtidaltensor_nonperiodic_diff(1);
#endif

    }
  if(i == 1)
    endrun(68686);
#ifdef DISTORTIONTENSORPS
//  check_tidaltensor_nonperiodic(10000);
#endif
#endif
#else
  i = pmforce(0);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
  pmtidaltensor_nonperiodic_fourier(0, 0);
  pmtidaltensor_nonperiodic_fourier(0, 1);
  pmtidaltensor_nonperiodic_fourier(0, 2);
  pmtidaltensor_nonperiodic_fourier(0, 3);
  pmtidaltensor_nonperiodic_fourier(0, 4);
  pmtidaltensor_nonperiodic_fourier(0, 5);
*/
/* FINITE DIFFERENCES based */
  pmtidaltensor_nonperiodic_diff(0);
#endif

  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      /* try again */
      i = pmforce(0);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
      pmtidaltensor_nonperiodic_fourier(0, 0);
      pmtidaltensor_nonperiodic_fourier(0, 1);
      pmtidaltensor_nonperiodic_fourier(0, 2);
      pmtidaltensor_nonperiodic_fourier(0, 3);
      pmtidaltensor_nonperiodic_fourier(0, 4);
      pmtidaltensor_nonperiodic_fourier(0, 5);
*/
/* FINITE DIFFERENCES based */
      pmtidaltensor_nonperiodic_diff(0);
#endif
    }
  if(i == 1)
    endrun(68687);
#ifdef DISTORTIONTENSORPS
//  check_tidaltensor_nonperiodic(10000);
#endif
#ifdef PLACEHIGHRESREGION
  i = pmforce(1);
#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
  pmtidaltensor_nonperiodic_fourier(1, 0);
  pmtidaltensor_nonperiodic_fourier(1, 1);
  pmtidaltensor_nonperiodic_fourier(1, 2);
  pmtidaltensor_nonperiodic_fourier(1, 3);
  pmtidaltensor_nonperiodic_fourier(1, 4);
  pmtidaltensor_nonperiodic_fourier(1, 5);
*/
/* FINITE DIFFERENCES based */
  pmtidaltensor_nonperiodic_diff(1);
#endif
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(i = 0; i < NumPart; i++)
	P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

#ifdef DISTORTIONTENSORPS
      P[i].tidal_tensorpsPM[0][0] = 0;
      P[i].tidal_tensorpsPM[0][1] = 0;
      P[i].tidal_tensorpsPM[0][2] = 0;
      P[i].tidal_tensorpsPM[1][0] = 0;
      P[i].tidal_tensorpsPM[1][1] = 0;
      P[i].tidal_tensorpsPM[1][2] = 0;
      P[i].tidal_tensorpsPM[2][0] = 0;
      P[i].tidal_tensorpsPM[2][1] = 0;
      P[i].tidal_tensorpsPM[2][2] = 0;
#endif
      i = pmforce(0) + pmforce(1);

#ifdef DISTORTIONTENSORPS
/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
      pmtidaltensor_nonperiodic_fourier(0, 0);
      pmtidaltensor_nonperiodic_fourier(0, 1);
      pmtidaltensor_nonperiodic_fourier(0, 2);
      pmtidaltensor_nonperiodic_fourier(0, 3);
      pmtidaltensor_nonperiodic_fourier(0, 4);
      pmtidaltensor_nonperiodic_fourier(0, 5);
*/
/* FINITE DIFFERENCES based */
      pmtidaltensor_nonperiodic_diff(0);

/* choose what kind of tidal field calculation you want */
/* FOURIER based */
/*
      pmtidaltensor_nonperiodic_fourier(1, 0);
      pmtidaltensor_nonperiodic_fourier(1, 1);
      pmtidaltensor_nonperiodic_fourier(1, 2);
      pmtidaltensor_nonperiodic_fourier(1, 3);
      pmtidaltensor_nonperiodic_fourier(1, 4);
      pmtidaltensor_nonperiodic_fourier(1, 5);
*/
/* FINITE DIFFERENCES based */
      pmtidaltensor_nonperiodic_diff(1);
#endif
    }
  if(i != 0)
    endrun(68688);
#ifdef DISTORTIONTENSORPS
//  check_tidaltensor_periodic(10000);
//  check_tidaltensor_nonperiodic(10000);
#endif
#endif
#endif


#ifndef PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0;

#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }


  /* Finally, the following factor allows a computation of cosmological simulation 
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#endif

}


#endif
