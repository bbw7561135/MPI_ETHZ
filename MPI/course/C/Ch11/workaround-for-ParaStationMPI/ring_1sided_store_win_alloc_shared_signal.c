/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the High Performance           *
 * Computing Centre Stuttgart (HLRS).                           *
 * The examples are based on the examples in the MPI course of  *
 * the Edinburgh Parallel Computing Centre (EPCC).              *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * HLRS and EPCC take no responsibility for the use of the      *
 * enclosed teaching material.                                  *
 *                                                              *
 * Authors: Joel Malard, Alan Simpson,            (EPCC)        *
 *          Rolf Rabenseifner, Traugott Streicher (HLRS)        *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                * 
 *                                                              *  
 * Purpose: A program to try out one-sided communication        *
 *          with window=rcv_buf and MPI_PUT to put              *
 *          local snd_buf value into remote window (rcv_buf).   *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

/* Courtesy to private email from Jed Brown, August 27, 2014 */
int volatile_load(int *addr)
{ return *(volatile int *)addr;
}

void volatile_store(int *addr, int val)
{ *((volatile int *)addr) = val;
}
  

int main (int argc, char *argv[])
{
  int my_rank, size;
  int snd_buf;
  int *rcv_buf_ptr, *signal_A_ptr, *signal_B_ptr;
  int right, left;
  int sum, i;

  MPI_Win  win, win_signal_A, win_signal_B;
  MPI_Comm comm_sm; int size_comm_sm;


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;
/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
/* right = my_rank + 1;          */
/* if (right == size) right = 0; */
/* left = my_rank - 1;           */
/* if (left == -1) left = size-1;*/

  /* Create the window. */

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
  MPI_Comm_size(comm_sm, &size_comm_sm); 
  if (size_comm_sm != size) { printf("Not on one shared memory node \n"); MPI_Abort(MPI_COMM_WORLD, 0); }

  /* ParaStation MPI may no allow MPI_Win_allocate_shared on MPI_COMM_WORLD. Workaround: Substitute MPI_COMM_WORLD by comm_sm (on 3 lines): */ 
  MPI_Win_allocate_shared((MPI_Aint)(1*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&rcv_buf_ptr, &win);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win);

  MPI_Win_allocate_shared((MPI_Aint)(1*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_A_ptr, &win_signal_A);
  *signal_A_ptr = 0;
  MPI_Win_fence(0,win_signal_A);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_A);

  MPI_Win_allocate_shared((MPI_Aint)(1*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_B_ptr, &win_signal_B);
  *signal_B_ptr = 0;
  MPI_Win_fence(0,win_signal_B);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_B);

  MPI_Win_sync(win); /* called by the origin process */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_sync(win); /* called by the target process */

  sum = 0;
  snd_buf = my_rank;

  for( i = 0; i < size; i++) 
  {
/*  ... The local Win_syncs are needed to sync the processor and real memory. */
/*  ... The following pair of syncs is needed that the read-write-rule is fulfilled. */
/*  rcv_buf is now exposed, i.e., the reading of previous value is finished */
    MPI_Win_sync(win); /* called by the target process */

/*  ... posting to left that rcv_buf is exposed to left, i.e.,
        the left process is now allowed to store data into the local rcv_buf  */
    volatile_store(signal_A_ptr+(left-my_rank), 1);
    while (volatile_load(signal_A_ptr)==0) /*IDLE*/;
    volatile_store(signal_A_ptr, 0);
/*  ... the two syncs for the both signals  can be removed, because
        the following  MPI_Win_sync(win) will sync all three variables */
/*  MPI_Win_sync(win_signal_A); */
/*  MPI_Win_sync(win_signal_B); */

    MPI_Win_sync(win); /* called by the origin process */
/*  rcv_buf is now writable by the origin process */     


    /* MPI_Put(&snd_buf, 1, MPI_INT, right, (MPI_Aint) 0, 1, MPI_INT, win); */
    /*   ... is substited by (with offset "right-my_rank" to store into right neigbor's rcv_buf): */
    *(rcv_buf_ptr+(right-my_rank)) = snd_buf; /* the origin process writes into the target rcv_buf */

/*  ... The following pair of syncs is needed that the write-read-rule is fulfilled. */
/*  writing of rcv_buf from remote is finished */
    MPI_Win_sync(win); /* called by the origin process */

/*  ... The following communication synchronizes the processors in the way */
/*      that the origin processor has finished the store */
/*      before the target processor starts to load the data.   */
/*  ... posting to right that rcv_buf was stored from left */
    volatile_store(signal_B_ptr+(right-my_rank), 1);
    while (volatile_load(signal_B_ptr)==0) /*IDLE*/;
    volatile_store(signal_B_ptr, 0);
/*  ... the two syncs for the both signals  can be removed, because
        the following  MPI_Win_sync(win) will sync all three variables */
/*  MPI_Win_sync(win_signal_B); */
/*  MPI_Win_sync(win_signal_A); */

    MPI_Win_sync(win); /* called by the target process */
/*  rcv_buf is now locally (i.e., by the origin process) readable */ 
    
    snd_buf = *rcv_buf_ptr;
    sum += *rcv_buf_ptr;
  }

  printf ("PE%i:\tSum = %i\n", my_rank, sum);

  MPI_Win_unlock_all(win);
  MPI_Win_unlock_all(win_signal_A);
  MPI_Win_unlock_all(win_signal_B);
  MPI_Win_free(&win);

  MPI_Finalize();
}
