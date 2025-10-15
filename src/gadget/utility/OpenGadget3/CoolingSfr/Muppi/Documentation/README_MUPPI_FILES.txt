These are the files that had to be changed to port MUPPI in OpenGadget

CoolingSfr/Muppi:
	muppi_communications.c
	muppi.h
	muppi_kicks.h
	muppi_normalize.c
	muppi_proto.h
	sfr_muppi.c

CoolingSfr/Sfr_LT:
	lt_sfr.c
	lt_sfr_light.c
	lt_sn.c
	lt_sn_light.c
	lt_sn.h
	lt_sn_init.c

CodeBase:
	allvars.c
	allvars.h
	begrun.c
        init.c
	proto.h
	run.c

Hydro:
	density.c
	hydra.c
        smooth_simple.c

IO:
	io.c
        read_ic.c
	
Integrator:
	kicks.c

System:
	tags.h



=======================================================================
This is the first working version of the porting of Muppi in 
OpenGadget, I'll name it OpenMuppi 0.1
26/04/2021

Verified OpenMuppi behaviour, Dec 2021, OpenMuppi 1.0

Ported Milena Valentini switches, OpenMuppi 1.1
10/2/2022
=======================================================================


=======================================================================
in this directory you will also find a Config.sh and a paramfile.txt
they are MUPPI-ready and intended to be run on AqC6 as IC
(not accluded here because it is relatively heavy)
You'll also find a paranoid-saved Makefile


With the porting of Milena Valentini switches, we added a file named
Config-MV-full.sh that includes them
=======================================================================


=======================================================================
TODO LIST
=======================================================================

1. Very urgent
   - Operational test of this porting   DONE
   - Cleaning up of remaining zombie code/comments

2. Urgent
   - Clean up the code in OpenGadget style
     (i.e. panici exits, verbose printouts, paramfile switches, etc)
   - Decision on the normalization round
     (standalone or in hydra?)
   - If we want to re-include in hydra, implement it
   - Detach MUPPI variables from SphP structure

3. Immediately following
   - Implement Milena's AGN feedback (as a switch) DONE
   - Merge with GL/Cinthia's dust 

4. Future
   - Re-implement all of our kinetic FB schemes
   - ...other things? suggestions?




