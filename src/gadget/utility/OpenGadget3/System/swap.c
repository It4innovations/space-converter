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
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include "../CodeBase/utilities.h"
#include "swap.h"

double SwapDouble(double Val)
{
  double nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

float SwapFloat(float Val)
{
  float nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

int SwapInt(int Val)
{
  int nVal;
  int i;
  const char *readFrom = (const char *) &Val;
  char *writeTo = ((char *) &nVal) + sizeof(nVal);

  for(i = 0; i < sizeof(Val); ++i)
    {
      *(--writeTo) = *(readFrom++);
    }
  return nVal;
}

int CheckSwap(char *fname, int *swap)
{
  FILE *fd;
  off_t fsize, fpos;
  int blocksize, blockend;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't open file `%s'.\n", fname);
      return -1;
    }

  fseeko(fd, 0, SEEK_END);
  fsize = ftello(fd);

  *swap = 0;
  fpos = 0;
  fseeko(fd, 0, SEEK_SET);
  safe_fread(&blocksize, sizeof(int), 1, fd);
  while(!feof(fd))
    {
      if(fpos + blocksize + 4 > fsize)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4 + blocksize;
      fseeko(fd, fpos, SEEK_SET);
      safe_fread(&blockend, sizeof(int), 1, fd);
      if(blocksize != blockend)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4;
      if(!fread(&blocksize, sizeof(int), 1, fd))
	break;
    }

  if(*swap == 0)
    {
      fclose(fd);
      return 0;
    }

  fpos = 0;
  fseeko(fd, 0, SEEK_SET);
  safe_fread(&blocksize, sizeof(int), 1, fd);
  while(!feof(fd))
    {
      blocksize = SwapInt(blocksize);
      if(fpos + blocksize + 4 > fsize)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4 + blocksize;
      fseeko(fd, fpos, SEEK_SET);
      safe_fread(&blockend, sizeof(int), 1, fd);
      blockend = SwapInt(blockend);
      if(blocksize != blockend)
	{
	  *swap += 1;
	  break;
	}
      fpos += 4;
      if(!fread(&blocksize, sizeof(int), 1, fd))
	break;
    }

  fclose(fd);
  return 0;
}
