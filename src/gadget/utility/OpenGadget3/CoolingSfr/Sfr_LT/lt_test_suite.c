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

/*-------------------------------- */
/*
  to be added in main.c right before the call to run()
*/

/* #ifdef LT_STELLAREVOLUTION */
/*   int i; */
/*   if(All.TestSuite) */
/*     { */
/*       i = get_random_number(100) * NTask; */
/*       if(i == ThisTask) */
/*         test_suite(); */

/*       MPI_Barrier(MYMPI_COMM_WORLD); */

/*       endrun(7777); */
/*     } */
/* #endif */


#if defined(LT_METAL_COOLING)
void test_metal_cooling()
{
  int i, j;

  double logT, logZ, rate;

  FILE *testfile;

  testfile = fopen("MetalCoolingRates.dat", "w");

  fprintf(testfile, "# metallicity (log solar), temperature (log K), rate (erg/sec*cm^3)\n#\n");
  for(i = 1; i < ZBins; i++)
    {
      logZ = (CoolZvalue[i - 1] + CoolZvalue[i]) / 2;
      for(j = 1; j < TBins; j++)
	{
	  logT = (CoolTvalue[j - 1] + CoolTvalue[j]) / 2;
	  rate = GetMetalLambda(logT, logZ);
	  fprintf(testfile, "%4.3e %5.4e %g\n", logZ, logT, rate);
	}
      fprintf(testfile, "#\n#\n");
    }
  fclose(testfile);
  return;
}
#endif



int test_suite()
{
  int test_result = 0, singletest_result;

  int i;

  double test1, test2;

  FILE *tfile, *file;


  tfile = fopen("TESTS.txt", "w");

  fprintf(tfile,
	  "\n-----------------------------------------------------------\n"
	  "TESTS made by processor %d\n\n", ThisTask);


  /* ----------------------------------------------------------------------------
   * test 1
   * test the cooling: this produces an ascii file, MetalCoolingRates.dat that
   * contains the cooling functions for temperature and metallcity values that
   * are within the table boundaries but not exactly equal to the tabulated values.
   * look in test_metal_cooling() for details about the file format.
   */
#if defined(LT_METAL_COOLING)
  fprintf(tfile, "test metal cooling rates.. ");
  test_metal_cooling();
  fprintf(tfile, "done (tables in MetalCoolingRates.dat)\n");
#endif


  /* ----------------------------------------------------------------------------
   * test 2
   * test all the IMF-related functions
   */
  /* this test should be made outside the code */




  /* ----------------------------------------------------------------------------
   * test 3
   * test the chemical evolution calculations
   */

  fprintf(tfile, "test chemical evolution.. ");
  TestStellarEvolution( );
  fprintf(tfile, "done (tables in SE.dbg)\n");


  /* ----------------------------------------------------------------------------
   * test 4
   * test the chemical time-stepping
   */

  /* to be done */


  /* ----------------------------------------------------------------------------
   * test 5
   * test the lifetime funcion and its inverse
   */

  fprintf(tfile, "test the lifetime function for 100 ranodm values.. ");
  for(i = 0; i < 100; i++)
    {
      test1 = All.mean_lifetime + (All.sup_lifetime - All.inf_lifetime) * get_random_number(i * 10);
      test2 = lifetime(dying_mass(test1));

      if(test2 <= 0)
	{
	  printf(">> WARN :: strange result for lifetime(dying_mass(T)):: T = %8.6e, res = %8.6e\n",
		 test1, test2);
	  singletest_result = 1;
	}
      else
	{
	  if(fabs(test2 - test1) / test2 > 1e-3)
	    {
	      printf(">> WARN :: inaccurate lifetime / lifetime^-1 for T = %8.6e"
		     "           lifetime^-1(T)           = %8.6e\n"
		     "           lifetime(lifetime^-1(T)) = %8.6e\n", test1, dying_mass(test1), test2);
	      singletest_result = 1;
	    }
	}
    }


  file = fopen("lifetime.txt", "w");
  test1 = log10(All.sup_lifetime / All.inf_lifetime) / 100;
  fprintf(file, "# time (Gyr) , dying_mass(time), lifetime(dying_amss(time)), dm(time)/dt\n");
  for(i = 0; i < 100; i++)
    {
      test2 = All.inf_lifetime * pow(10, test1 * i);
      fprintf(file, "%8.6e\t%8.6e\t%8.6e\t%8.6e\n",
	      test2, dying_mass(test2), lifetime(dying_mass(test2)), dm_dt(dying_mass(test2), test2));
    }
  fclose(file);

  test_result += singletest_result;
  fprintf(tfile, "done\n");



  fclose(tfile);
  return test_result;
}



void TestStellarEvolution( )
{
#define Mass_N 30
#define Z_N 8

  double top_masses[Mass_N+1], bot_masses[Mass_N+1];
  double mymetals[LT_NMet], myenergy, numsn, Zstar;
  FILE *outfile[3];
  char *names[3] = {"SnIa","SnII","LIMS"};
  int   idx;
  
  struct particle_data     myP, *save_P;
  struct met_particle_data myMetP, *save_MetP;

  printf("Start testing..\n");

  /*  as first, we integrate SnII, SnIa and AGB on
   *  a mass grid; the results will be compared
   *  off-line with external routines to validete them
  */

  top_masses[0] = 50;
  bot_masses[0] = 35;
  top_masses[Mass_N - 1] = 1.0;
  bot_masses[Mass_N - 1] = 0.8;
      
  double max = log10(35);
  double min = log10(0.6);
  double delta = (max - min)/Mass_N;
      
  for( int i = 0; i <= Mass_N - 1; i++)
    {
      bot_masses[i] = pow(10, min + delta*i);
      top_masses[i] = bot_masses[i]*1.2;
    }

  save_P    = P;
  P         = &myP;
  save_MetP = MetP;
  MetP      = &myMetP;

  P[0].Type = 4;
  P[0].pt.MetID = 0;
  P[0].Mass = 1.0;

  MetP[0].PID = 0;
  MetP[0].iMass = 1;
  for( int i = 0; i < LT_NMet; i++ )
    MetP[0].Metals[i] = 0;
  MetP[0].Metals[Hel] = (1 - HYDROGEN_MASSFRAC);
  idx = (Hel > 0 ? 0 : 1);
  
  for ( int tt = 0; tt < 3; tt++ )
    {
      char buffer[100];
      sprintf( buffer, "SE_%s.dbg", names[tt] );

      if((outfile[tt] = fopen(buffer, "w")) == 0x0)
	{
	  printf("it has been impossible to open file %s..\n", buffer);
	  fflush(stdout);
	  return;
	}
    }


  delta = log10(0.2 / 2e-5) / (Z_N - 1);

  if(UseSnII)
    {
      fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
      fprintf(FdSnInit, "effective yields on a Z grid - SnII from %g to %g\n",
	      dying_mass(All.mean_lifetime), dying_mass(SFs[0].ShortLiv_TimeTh));
      for( int k = 0; k < Z_N; k++ )
	{
	  Zstar = 2e-5 * pow(10, delta * k);
	  MetP[0].Metals[idx] = Zstar / (1 + Zstar) * 0.76;
	  
	  numsn = 0;
	  get_SnII_product(0, 0, &mymetals[0], &myenergy, dying_mass(All.mean_lifetime), dying_mass(SFs[0].ShortLiv_TimeTh), &numsn);
	  fprintf(FdSnInit, "%10.8e %10.8e %10.8e :: %10.8e :: ",
		 Zstar, SFs[0].ShortLiv_TimeTh, All.mean_lifetime, numsn);
	  for ( int j = 0; j < LT_NMet; j++ )
	    fprintf(FdSnInit, "%10.8e ", mymetals[j]);
	  fprintf(FdSnInit, "\n"); fflush(FdSnInit);
	}      
    }

  if(UseSnIa)
    {
      double sup = lifetime(All.MBms)*0.999;
      fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
      fprintf(FdSnInit, "effective yields on a Z grid - SnIa from %g %g\n",
	      All.mean_lifetime, sup);
      
      get_SnIa_product(0, 0, &mymetals[0], &myenergy, All.mean_lifetime, sup);
      fprintf(FdSnInit, "%10.8e %10.8e %10.8e :: - :: ",
	     Zstar, All.mean_lifetime, sup);
      for ( int j = 0; j < LT_NMet; j++ )
	fprintf(FdSnInit, "%10.8e ", mymetals[j]);
      fprintf(FdSnInit, "\n"); fflush(FdSnInit);
    }

  if(UseAGB)
    {
      fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
      fprintf(FdSnInit, "effective yields on a Z grid - AGB from %g %g\n",
	      All.MBms, All.Mup);
      for( int k = 0; k < Z_N; k++ )
	{
	  Zstar = 2e-5 * pow(10, delta * k);
	  MetP[0].Metals[idx] = Zstar / (1 + Zstar) * 0.76;
	  
	  numsn = 0;
	  double sup = lifetime(All.MBms)*0.999;
	  get_AGB_product(0, 0, &mymetals[0], All.MBms, All.Mup, &numsn);
	  fprintf(FdSnInit, "%10.8e %10.8e %10.8e %10.8e :: - :: ",
		 Zstar, All.MBms, All.Mup, numsn);
	  for ( int j = 0; j < LT_NMet; j++ )
	    fprintf(FdSnInit, "%10.8e ", mymetals[j]);
	  fprintf(FdSnInit, "\n"); fflush(FdSnInit);
	}
    }

  printf("\n");

  
  printf("SSP evolution on a time grid\n");
  
  for( int k = 0; k < Z_N; k++ )
    {
      fflush(stdout);
      Zstar = 2e-5 * pow(10, delta * k);
      MetP[0].Metals[idx] = Zstar / (1 + Zstar) * 0.76;


      for( int i = 0; i <= Mass_N - 1; i++ )
	{
	  double tstart = lifetime(top_masses[i]);
	  double tend   = lifetime(bot_masses[i]);
	  double deltat = tend - tstart;
	  
	  if(UseSnII)
	    {
	      get_SnII_product(0, 0, &mymetals[0], &myenergy, bot_masses[i], top_masses[i], &numsn);
	      fprintf(outfile[1], "%10.8e %10.8e %10.8e %10.8e %10.8e :: %10.8e ",
		      Zstar, tstart, tend, bot_masses[i], top_masses[i],
		      numsn);
	      for ( int j = 0; j < LT_NMet; j++ )
		fprintf(outfile[1], "%10.8e ", mymetals[j]);
	      fprintf(outfile[1], "\n");
	    }


	  if(UseSnIa)
	    {
	      gsl_function myF;
	      double abserr;
	      myF.function = &nRSnIa;
	      myF.params = &IMFs[0];

	      if( tend > All.mean_lifetime )
		{
		  double mytstart = (tstart > All.mean_lifetime ? tstart : All.mean_lifetime);
		  double mydeltat = tend - mytstart;
		  gsl_integration_qag(&myF, mytstart, tend,
				      qag_ABS_ERR, qag_ABS_ERR, 1000, qag_INT_KEY, w, &numsn, &abserr);
		  numsn *= MetP[0].iMass;
		  get_SnIa_product(0, 0, &mymetals[0], &myenergy, mytstart, mydeltat );
		}
	      else
		{
		  numsn = 0;
		  memset( mymetals, 0, sizeof(double)*LT_NMet);
		}
	      fprintf(outfile[0], "%10.8e %10.8e %10.8e %10.8e %10.8e :: %10.8e ", Zstar,
		      tstart, tend, bot_masses[i], top_masses[i], numsn);
	      for( int j = 0; j < LT_NMet; j++ )
		fprintf(outfile[0], "%10.8e ", mymetals[j]);
	      fprintf(outfile[0], "\n");
	    }

	  if(UseAGB)
	    {
	      get_AGB_product(0, 0, &mymetals[0], bot_masses[i], top_masses[i], &numsn);
		
	      fprintf(outfile[2], "%10.8e %10.8e %10.8e %10.8e %10.8e :: %10.8e ",
		      Zstar, tstart, tend, bot_masses[i], top_masses[i],
		      numsn);
	      for( int j = 0; j < LT_NMet; j++ )
		fprintf(outfile[2], "%10.8e ", mymetals[j]);
	      fprintf(outfile[2], "\n");
	    }
	}
      for ( int tt = 0; tt < 3; tt++ )
	fprintf(outfile[tt], "\n\n");
    }
  for ( int tt = 0; tt < 3; tt++ )
    fclose(outfile[tt]);



  printf("impact of round-off on SnII\n");
  /*
   * Now let's check the impact of round-off error
   * in floating-point summation
   */

  MetP[0].Metals[idx] = 0;

  float  SnII_sum[LT_NMet] = {0};
  double SnII_dsum[LT_NMet] = {0};
  double SnII_all[LT_NMet] = {0};

  fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
  fprintf(FdSnInit, "effect of round-off on SnII from %g to %g\n",
	  dying_mass(All.mean_lifetime), dying_mass(SNtimesteps[0][0][1])); fflush(stdout);
  
  get_SnII_product(0, 0, &SnII_all[0], &myenergy,
		    dying_mass(All.mean_lifetime), dying_mass(SNtimesteps[0][0][1]), &numsn);
  for ( int i = 1; i < ShortLiv_Nsteps[0]; i++ )
    {
      double metals[LT_NMet] = {0};
      double sup_mass = dying_mass(SNtimesteps[0][0][i]);
      double inf_mass = dying_mass(SNtimesteps[0][0][i]+SNtimesteps[0][1][i]);
      get_SnII_product(0, 0, &metals[0], &myenergy, inf_mass, sup_mass, &numsn);
      for( int j = 0; j < LT_NMet; j++ )
	SnII_sum[j] += metals[j], SnII_dsum[j] += metals[j];
    }


  for( int j = 0; j < LT_NMet; j++ )
    fprintf(FdSnInit, "\t element %d : %g %g %g (%g)\n", j,
	    SnII_all[j], SnII_dsum[j], SnII_sum[j], SnII_all[j] > 0 ? (SnII_sum[j] - SnII_all[j])/SnII_all[j] : 0);

 
  int last = Nsteps[0]-2;
  printf("impact of round-off on SnIa\n");
  float  SnIa_sum[LT_NMet] = {0};
  double SnIa_dsum[LT_NMet] = {0};
  double SnIa_all[LT_NMet] = {0};  

  fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
  fprintf(FdSnInit, "effect of round-off on SnIa from %g to %g\n", SNtimesteps[0][0][ShortLiv_Nsteps[0]], SNtimesteps[0][0][last]+SNtimesteps[0][1][last]); fflush(stdout);

  get_SnIa_product(0, 0, &SnIa_all[0], &myenergy,
		    SNtimesteps[0][0][ShortLiv_Nsteps[0]],
		    SNtimesteps[0][0][last]+SNtimesteps[0][1][last]);
  for ( int i = ShortLiv_Nsteps[0]; i <= last; i++ )
    {
      double metals[LT_NMet] = {0};
      double inf_time = SNtimesteps[0][0][i];
      double sup_time = SNtimesteps[0][0][i]+SNtimesteps[0][1][i];
      get_SnIa_product(0, 0, &metals[0], &myenergy, inf_time, sup_time-inf_time);
      for( int j = 0; j < LT_NMet; j++ )
	SnIa_sum[j] += metals[j], SnIa_dsum[j] += metals[j];
      printf("\n");
    }

  
  for( int j = 0; j < LT_NMet; j++ )
    fprintf(FdSnInit, "\t element %d : %g %g %g (%g)\n", j,
	    SnIa_all[j], SnIa_dsum[j], SnIa_sum[j], SnIa_all[j] > 0 ? (SnIa_sum[j] - SnIa_all[j])/SnIa_all[j] : 0);
  fprintf(FdSnInit, "\n");


  printf("impact of round-off on AGB\n");
  float  AGB_sum[LT_NMet] = {0};
  double AGB_dsum[LT_NMet] = {0};
  double AGB_all[LT_NMet] = {0};  

  fprintf(FdSnInit, "\n----------------------------------------------------------\n\n");
  fprintf(FdSnInit, "effect of round-off on AGB from %g to %g\n",
	  SNtimesteps[0][0][ShortLiv_Nsteps[0]], SNtimesteps[0][0][last]+SNtimesteps[0][1][last]); fflush(stdout);

  get_AGB_product(0, 0, &AGB_all[0],
		    SNtimesteps[0][0][ShortLiv_Nsteps[0]],
		   SNtimesteps[0][0][last]+SNtimesteps[0][1][last], &numsn);
  for ( int i = ShortLiv_Nsteps[0]; i <= last; i++ )
    {
      double metals[LT_NMet] = {0};
      double inf_time = SNtimesteps[0][0][i];
      double sup_time = SNtimesteps[0][0][i]+SNtimesteps[0][1][i];
      get_AGB_product(0, 0, &metals[0], inf_time, sup_time, &numsn);
      for( int j = 0; j < LT_NMet; j++ )
	AGB_sum[j] += metals[j], AGB_dsum[j] += metals[j];
    }
  
  for( int j = 0; j < LT_NMet; j++ )
    fprintf(FdSnInit, "\t element %d : %g %g %g (%g)\n", j,
	    AGB_all[j], AGB_dsum[j], AGB_sum[j], AGB_all[j] > 0 ? (AGB_sum[j] - AGB_all[j])/AGB_all[j] : 0);
  fprintf(FdSnInit, "\n");
    
  fflush( FdSnInit );

  P    = save_P;
  MetP = save_MetP;

  return;
}


void TestStarsEvolution( int N, int snap )
{
  FILE *myfile;
  char fname[200];
  long double array[LT_NMetP] = {0};


  sprintf( fname, "stars_evolution.%d.%s_%d", ThisTask, (snap < 0 ? "tstep" : "snap"), (snap < 0 ? All.NumCurrentTiStep : snap));
  myfile = fopen( fname, "w");
  if ( myfile == NULL )
    return;

  fprintf(myfile, "%g %g %g\n", All.Time, All.Time_Age, get_age(All.Time));

  N = ( N  < 0 ? N_stars : N );
  
  for ( int i = 0; i < N; i++ )
    {
      int    mySFi, myIMFi, Yset;
      
      /* find the IMF associated to this particle */
      get_SF_index(MetP[i].PID, &mySFi, &myIMFi);
      
      /* find the yield set associated to this IMF */
      Yset = IMFs[myIMFi].YSet;
      
      double mylifetime;
      mylifetime = get_age(MetP[i].StellarAge) - All.Time_Age;
      
      if ( mylifetime <= SFs[mySFi].ShortLiv_TimeTh )
	continue;
      
      
      double metals[3][LT_NMet] = {0}, all_metals[LT_NMet];
      double inf[3] = {0}, sup[3] = {0};
      
      for( int s = 0; s < 3; s++ )
	for(int j = 0; j < LT_NMet; j++)
	  metals[s][j] = 0;


      /* type II */
      inf[0] = (mylifetime > All.mean_lifetime ? All.Mup : dying_mass(mylifetime));
      sup[0] = dying_mass(SFs[mySFi].ShortLiv_TimeTh);
      double energy, numsn;
      double mass = get_SnII_product(MetP[i].PID, Yset, &metals[0][0], &energy, inf[0], sup[0], &numsn);
      
      for(int j = 0; j < LT_NMet; j++)
	array[j] += metals[0][j];

      /* Ia and AGB */
      if ( mylifetime > All.mean_lifetime )
	{
	  inf[1] = All.mean_lifetime;
	  sup[1] = mylifetime;

	  double energy;
	  double mass = get_SnIa_product(MetP[i].PID, Yset, &metals[1][0], &energy, inf[1], mylifetime - All.mean_lifetime);

	  for(int j = 0; j < LT_NMet; j++)
            array[j] += metals[1][j]; 

	  inf[2] = dying_mass(mylifetime);
	  sup[2] = All.Mup;

	  double numsn;
	  mass = get_AGB_product(MetP[i].PID, Yset, &metals[2][0], inf[2], sup[2], &numsn);

	  for(int j = 0; j < LT_NMet; j++) 
	    array[j] += metals[2][j];	  
	}

      double Zstar = get_metallicity(i, -1);
      fprintf( myfile, "%llu %lg %lg %lg %lg %lg - ", (unsigned long long)P[MetP[i].PID].ID, MetP[i].iMass, Zstar, MetP[i].StellarAge, get_age(MetP[i].StellarAge), mylifetime);
      double sumZ = 0;
      for(int j = 0; j < LT_NMetP; j++){
	sumZ += MetP[i].Metals[j];
	fprintf(myfile, "%lg ", MetP[i].Metals[j] );}
      fprintf(myfile, "%lg ", MetP[i].iMass - sumZ);
      fprintf(myfile, "- ");
      for ( int s = 0; s < 3; s++ )
	fprintf(myfile, "%lg %lg ", inf[s], sup[s] );
      fprintf(myfile, "- ");
      for ( int s = 0; s < 3; s++ ) {
	for(int j = 0; j < LT_NMet; j++)
	  fprintf(myfile, "%lg ", metals[s][j] );
	fprintf(myfile, ": ");}
      fprintf(myfile, "\n");
    }

  fclose(myfile);

  MPI_Barrier(MPI_COMM_WORLD);

  long double total_array[NTask * LT_NMetP];
  memset( total_array, 0, NTask * LT_NMetP * sizeof(long double) );
  MPI_Gather(array, LT_NMetP*sizeof(long double), MPI_BYTE, total_array, LT_NMetP*sizeof(long double), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if ( ThisTask ==  0)
    {
      for( int j = 0; j < LT_NMetP; j++ )
	{
	  for( int i = 1; i < NTask; i++ )
	    total_array[j] += total_array[i*LT_NMetP + j];
	  array[j] = total_array[j];
	}

      printf("CHECKSUM stellarevolution %d %g : ", All.NumCurrentTiStep, All.Time);
      for( int j = 0; j < LT_NMetP; j++ )
	printf("%Lg ", total_array[j] );
      printf("\n");
    }

  MPI_Barrier(MYMPI_COMM_WORLD);

  return;
}


void SimulateStarsEvolution( )
{
  #define Ns 3
  struct particle_data     myP, *save_P;
  struct met_particle_data myMetP, *save_MetP;

  /* float imass[Ns] = {0.000156539, 0.000156504}; */
  /* double imetals[Ns][LT_NMetP] = {{3.75775e-05, 8.93066e-10, 5.14095e-11,1.03086e-08,5.43524e-12,1.88759e-09,6.82222e-10,4.13497e-10,9.65949e-10,3.57118e-10,5.02599e-12,1.55895e-11,6.52333e-11,1.15284e-11,5.2115e-13}, */
  /* 				  {3.75658e-05, 4.66233e-10,2.27601e-11, 4.55991e-09, 7.49774e-12, 9.18248e-10, 3.85286e-10, 1.83483e-10, 4.55286e-10, 2.35313e-10, 3.78572e-12, 1.12141e-11, 2.84799e-11, 7.4496e-12,1.27325e-13} }; */
  /* double inf[Ns][3] = {{8,0.0286133,4.38966},{8,0.0286133,7.70868}}; */
  /* double sup[Ns][3] = {{40,0.0949062,8},{40,0.0304327,8}}; */
  

  // 25.snap_6
  float imass[Ns] = {0.000156632, 0.000156472};
  double imetals[Ns][LT_NMetP] = { {3.76097e-05,1.68631e-09,7.13953e-11,1.36783e-08,1.77732e-11,2.70834e-09,1.34304e-09,5.75542e-10,1.51341e-09,9.01054e-10,1.10545e-11, 4.02939e-11, 8.80354e-11, 2.70998e-11, 0},
				   {3.75541e-05, 2.29196e-10, 3.35326e-11, 1.10141e-09, 3.6045e-12, 1.13074e-10, 1.34142e-10, 2.43486e-10, 4.75325e-10,1.90267e-09,6.90713e-13,6.09802e-12,4.29761e-11,2.71811e-10,3.20911e-11 } };
				  
  double inf[Ns][3] = {{8,0.0286133,4.2534},{8,0.0286133, 5.57507}};
  double sup[Ns][3] = {{40,0.101858,8},{40,0.0564619,8}};

  int   idx = !Hel;
  
  save_P    = P;
  P         = &myP;
  save_MetP = MetP;
  MetP      = &myMetP;

  P[0].Type = 4;
  P[0].pt.MetID = 0;
  MetP[0].PID = 0;

  for( int n = 0; n < Ns; n++ )
    {
      MetP[0].iMass = imass[n];
      for ( int j = 0; j < LT_NMetP; j++ )
	MetP[0].Metals[j] = imetals[n][j];

      double metals[3][LT_NMet];
      double energy, numsn;
      double mass;
      
      mass = get_SnII_product(0, 0, &metals[0][0], &energy, inf[n][0], sup[n][0], &numsn);
      mass = get_SnIa_product(0, 0, &metals[1][0], &energy, inf[n][1], sup[n][1]-All.mean_lifetime);
      mass = get_AGB_product(0, 0, &metals[2][0], inf[n][2], sup[n][2], &numsn);
      
    }  

  MetP = save_MetP;
  P    = save_P;

  return;
}
  
void get_metals_checksum( int mode, long double *checksum )
{
  long double local_checksum[LT_NMetP];
  memset( local_checksum, 0, LT_NMetP * sizeof(long double) );
  
  for( int i = 0; i < N_gas; i++ )
    for( int j = 0; j < LT_NMetP; j++ )
      local_checksum[j] += ( j!= Hyd ? SphP[i].Metals[j] : 0);

  if ( mode )
    {
      for ( int i = N_gas; i < NumPart; i++ )
	switch( P[i].Type )
	  {
	  case 4: { int idx = P[i].pt.MetID;
	      for( int j = 0; j < LT_NMetP; j++ )
		local_checksum[j] += ( j!=Hyd ? MetP[idx].Metals[j] : 0); }
	    break;
	  /* case 5: { int idx = P[i].BHID; */
	  /*     for( int j = 0; j < LT_NMetP; j++ ) */
	  /* 	local_checksum[j] += BhP[idx].Metals[j]; } */
	  /*   break; */
	  default: break;
	  }
    }

  long double total_checksum[NTask * LT_NMetP];
  memset( total_checksum, 0, NTask * LT_NMetP * sizeof(long double) );
  MPI_Gather(local_checksum, LT_NMetP*sizeof(long double), MPI_BYTE, total_checksum, LT_NMetP*sizeof(long double), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if ( ThisTask ==  0)
    for( int j = 0; j < LT_NMetP; j++ )
      {
	for( int i = 1; i < NTask; i++ )
	  total_checksum[j] += total_checksum[i*LT_NMetP + j];
	checksum[j] = total_checksum[j];
      }

  MPI_Barrier(MYMPI_COMM_WORLD);
  return;
}
