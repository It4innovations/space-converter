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
/** \file
    MPI utility functions.
*/

#include <mpi.h>
#include <string.h>
#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"
#include "string.h"


/** Calculates the recv_count, send_offset, and recv_offset arrays
    based on the send_count. Returns nimportbytes, the total number of
    bytes to be received. If an identical set of copies are to be
    sent to all tasks, set send_identical=1 and the send_offset will
    be zero for all tasks.

    All arrays should be allocated with NTask size. */
int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset,
			  int send_identical)
{
  // Exchange the send/receive counts
  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);

  int nimportbytes = 0;
  recv_offset[0] = 0;
  send_offset[0] = 0;
  int j;
  for(j = 0; j < NTask; j++)
    {
      nimportbytes += recv_count[j];

      if(j > 0)
	{
	  send_offset[j] = send_offset[j - 1] + (send_identical ? 0 : send_count[j - 1]);
	  recv_offset[j] = recv_offset[j - 1] + recv_count[j - 1];
	}
    }
  return nimportbytes;
}


/** Compare function used to sort an array of int pointers into order
    of the pointer targets. */
int intpointer_compare(const void *a, const void *b)
{
  if((**(int **) a) < (**(int **) b))
    return -1;

  if((**(int **) a) > (**(int **) b))
    return +1;

  return 0;
}


/** Sort an opaque array into increasing order of an int field, given
    by the specified offset. (This would typically be field indicating
    the task.) Returns a sorted copy of the data array, that needs to
    be myfreed.

    We do this by sorting an array of pointers to the task field, and
    then using this array to deduce the reordering of the data
    array. Unfortunately this means making a copy of the data, but
    this just replaces the copy after the mpi_exchange_buffers
    anyway.  */
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr)
{
  int i;
  int **perm;
  *data2ptr = (char *) mymalloc_movable(data2ptr, "data2", n_items * item_size);
  perm = (int **) mymalloc("perm", n_items * sizeof(*perm));

  for(i = 0; i < n_items; ++i)
    perm[i] = (int *) ((char *) data + i * item_size + field_offset);

  qsort(perm, n_items, sizeof(*perm), intpointer_compare);
  // reorder data into data2
  for(i = 0; i < n_items; ++i)
    {
      size_t orig_pos = ((char *) perm[i] - ((char *) data + field_offset)) / item_size;
      if(!(((char *) perm[i] - ((char *) data + field_offset)) % item_size == 0))
	PANIC("something wrong here!");
      memcpy((char *) *data2ptr + item_size * i, (char *) data + item_size * orig_pos, item_size);
    }

  myfree(perm);
}


/** This function distributes the members in an opaque structure to
    the tasks based on a task field given by a specified offset into
    the opaque struct. The task field must have int type. n_items is
    updated to the new size of data. max_n is the allocated size of
    the data array, and is updated if a realloc is necessary.  */
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size)
{
  int i;

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < *n_items; i++)	/* find Send_count for each task */
    {
      int task = *((int *) ((char *) data + i * item_size + task_offset));
      if(!(task >= 0 && task < NTask))
	PANIC("wrong task number!");
      Send_count[task] += item_size;
    }

  void *data2;
  sort_based_on_field(data, task_offset, *n_items, item_size, &data2);	/* sort data based on target task */

  int nimportbytes = mpi_calculate_offsets(Send_count, Send_offset, Recv_count, Recv_offset, 0);	/* calculate offsets of send and receive buffers */
  if(nimportbytes % item_size != 0)
    PANIC("should not happen!");
  int nimport = nimportbytes / item_size;

  if(*max_n < nimport)		/* realloc data if too small */
    {
      data = (char *) myrealloc_movable(data, nimportbytes);
      *max_n = nimport;
    }

  MPI_Alltoallv(data2, Send_count, Send_offset, MPI_BYTE, data, Recv_count, Recv_offset, MPI_BYTE, MYMPI_COMM_WORLD);	/* exchange data */

  myfree(data2);

  *n_items = nimport;
}
