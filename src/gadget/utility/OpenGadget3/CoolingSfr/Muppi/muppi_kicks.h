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


#include "../../CodeBase/allvars.h"
#include "../../CodeBase/proto.h"


/*! \file muppi_kicks.c
*  \brief Add the kinetic energy for MUPPI kinetic feedback, called my: Integrator/kicks.c
*
*/
#ifdef GM_MUPPI

static inline void muppi_kinetic(int i)
{
  MyAtLeastDouble a3inv, vkick, vrsr;
  MyAtLeastDouble dx_i, dy_i, dz_i, vrsr_i, dir_i_norm;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1.0;


  /* THERMAL energy */
  /* note, here flows have already been considered */
  if(P[i].Type == 0)
    {
      SphP[i].Entropy += SphP[i].E_rec /	/* energy from MP parts */
	/* energy to specific energy */
	(P[i].Mass - SphP[i].M_sf) *
	/* energy to specific energy */
	GAMMA_MINUS1 / pow((double) SphP[i].Density * a3inv, GAMMA_MINUS1);
      /* ... to entropy */
      
      SphP[i].E_rec = 0.0;


      /* MUPPI: */
      /* velocity kick due to KINETIC energy */
      /* 
	 ENERGY scheme
	 ALL energy (summed as a scalar) is given in the direction corresponding to the weighted average of vector sum
	 of directions j-i. The weight is the kinetic energy. */
      
      if(SphP[i].E_kin > 0.0)
	{
	  vkick = sqrt(2.0 * SphP[i].E_kin / P[i].Mass);
	  if(All.ComovingIntegrationOn)
	    vkick *= All.Time;
	  
	  dx_i = -SphP[i].GradDens[0];
	  dy_i = -SphP[i].GradDens[1];
	  dz_i = -SphP[i].GradDens[2];
	  
	  dir_i_norm = SphP[i].GradDens[0]*SphP[i].GradDens[0]
	    + SphP[i].GradDens[1]*SphP[i].GradDens[1]
	    + SphP[i].GradDens[2]*SphP[i].GradDens[2];    /* milena */
	  
	  vrsr_i = sqrt(dir_i_norm);
	  
	  SphP[i].xkin  =  (dx_i/vrsr_i);   /* milena */
	  SphP[i].ykin  =  (dy_i/vrsr_i);
	  SphP[i].zkin  =  (dz_i/vrsr_i);

#ifdef MV_GM_STELLAR_KIN_FB2_OUTPUT
	  if( vkick/All.Time > 5000.0) {
	    printf("   V in kicks.c  time: %f   vkick: %e part ID: %d  hsml: %f     Temp: %e\n",
		   All.Time, vkick/All.Time, P[i].ID, P[i].Hsml, SphP[i].Temperature); fflush(stdout);
	  }
	  
	  fprintf(FdWind,"%g   %f %f %f   %f %f %f   %e   %e %e %e   %d   %e %e   %e %e %e   %e\n",
		  All.Time, P[i].Pos[0], P[i].Pos[0], P[i].Pos[2],
		  P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], P[i].Mass,
		  SphP[i].Density, SphP[i].Pressure, SphP[i].Temperature, P[i].ID,
		  vkick/All.Time, SphP[i].E_kin,
		  SphP[i].xkin, SphP[i].ykin, SphP[i].zkin,
		  SphP[i].DelayTime); 
	  fflush(FdWind);	
#endif


	  /* do the kick! ENERGY */
	  P[i].Vel[0] += vkick * SphP[i].xkin;
	  P[i].Vel[1] += vkick * SphP[i].ykin;
	  P[i].Vel[2] += vkick * SphP[i].zkin;
	  SphP[i].VelPred[0] += vkick * SphP[i].xkin;
	  SphP[i].VelPred[1] += vkick * SphP[i].ykin;
	  SphP[i].VelPred[2] += vkick * SphP[i].zkin;
	  
	  
	  SphP[i].Ekin_rec[0] += vkick * SphP[i].xkin;
	  SphP[i].Ekin_rec[1] += vkick * SphP[i].ykin;
	  SphP[i].Ekin_rec[2] += vkick * SphP[i].zkin;
	  
	  
	}
      SphP[i].E_kin = 0.0;
      SphP[i].xkin = 0.0;
      SphP[i].ykin = 0.0;
      SphP[i].zkin = 0.0;
    }

}

#endif //closes GM_MUPPI
