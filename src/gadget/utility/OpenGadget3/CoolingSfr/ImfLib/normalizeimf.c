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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libimf_vars.h>

#define Chabrier_MC 0.079
#define Chabrier_Sigma2 (0.69*0.69)

double Chabrier_Exp(double);
double Chabrier_byNum(double, void *);
double Chabrier_byMass(double, void *);
double Chabrier_byEgy(double, void *);

gsl_integration_workspace *my_w;

int main(int argc, char **argv)
{
  char *IMFfilename;
  double inf, sup, A;
  int i, j;

  if(argc < 4)
    {
      printf("arguments: IMF(s) file name\n"
             "           inf mass limit\n"
             "           sup mass limit\n");
      return -1;
    }

  my_w = gsl_integration_workspace_alloc(limf_gsl_intspace_dim);

  IMFfilename = (char*)malloc(strlen(*(argv+1)) + 2);
  sprintf(IMFfilename, "%s", *(argv+1));

  inf = atof(*(argv+2));
  sup = atof(*(argv+3));

  initialize_externalIMFs(1);
  set_externalIMF(0, "Chabrier", &Chabrier_byMass, &Chabrier_byNum, 0x0);

  read_imfs(IMFfilename);
  IMFp = &IMFs[0];

  printf("\n");

  for(j = 0; j < IMFs_dim; j++)
    {
      printf("------------------------------------------------------\n");
      printf("IMF :: %s :: ", IMFs[j].name);

      A = IntegrateIMF_byMass(inf, sup, &IMFs[j], INC_BH);
      if( fabs(A-1)/A > IMFs[j].Mm/1000)
        {
          printf("renormalizing.. (%g) ", fabs(A-1)/A);
          A = Renormalize_IMF(&IMFs[j], inf, sup, INC_BH);
        }
      printf("\n\n");
      if(fabs(inf - IMFs[j].Mm)/inf > 1e-2)
        printf("  >> warning : inf normalization limit is different than inf mass of IMF\n");
      if(fabs(sup - IMFs[j].MU)/sup > 1e-2)
        printf("  >> warning : sup normalization limit is different than sup mass of IMF\n");

      for(i = IMFs[j].NSlopes-1; i>=0; i--)
        printf("\tmass range :: [%.2f -> %.2f]\n\t\tnormalization :: %.7f\n", IMFs[j].Slopes.masses[i], (i>0)?IMFs[j].Slopes.masses[i-1]:IMFs[j].MU, IMFs[j].A[i]);
      
      printf("\nnormalization now reads :: %g\n\n", 1.0/A);
    }
  
  for(j = 0; j < IMFs_dim; j++)
    print_IMF(j, IMFfilename);
  return 0;
}


double Chabrier_Exp(double arg)
{
  return exp(-(log10(arg) - log10(Chabrier_MC)) * (log10(arg) - log10(Chabrier_MC))) / (2 * Chabrier_Sigma2);
}
double Chabrier_byMass(double arg, void *param)
{

  if(arg > 1)
    return IMFp->A[0] * pow(arg, -IMFp->Slopes.slopes[0]);
  else
    return IMFp->A[1] * Chabrier_Exp(arg);
}

double Chabrier_byNum(double arg, void *param)
{

  if(arg > 1)
    return IMFp->A[0] * pow(arg, -(1 + IMFp->Slopes.slopes[0]));
  else
    return IMFp->A[1] * Chabrier_Exp(arg);
}

