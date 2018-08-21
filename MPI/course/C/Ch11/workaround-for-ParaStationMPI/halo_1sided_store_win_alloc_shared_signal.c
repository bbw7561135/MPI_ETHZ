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
 * Purpose: A program to meassure 1-dim halo communication      *
 *          in myrank -1 and +1 directions (left and right)     *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define number_of_messages 50
#define start_length 4
#define length_factor 8
#define max_length 8388608 /* ==> 2 x 32 MB per process */
#define number_package_sizes 8
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */

/* Distance of two signal variables (in number of ints): */
/* To prohibit cache-line false sharing, you may choose at least 32 with 128 byte cache lines */
#define signal_distance 32
/* Signal exchange with memory fences: yes=1, no=0 */
#define signal_fence 0

/* Courtesy to private email from Jed Brown, August 27, 2014 */
int volatile_load(MPI_Win win, int *addr)
{ 
#if(signal_fence==1)
  MPI_Win_sync(win);
#endif 
  return *(volatile int *)addr;
}

void volatile_store(MPI_Win win, int *addr, int val)
{ *((volatile int *)addr) = val;
#if(signal_fence==1)
  MPI_Win_sync(win);
#endif 
}


int main(int argc, char *argv[])
{
  int i, j, k, length, my_rank, left, right, size, test_value, mid;    
  double start, finish, transfer_time; 
  float snd_buf_left[max_length], snd_buf_right[max_length];
  float *rcv_buf_left, *rcv_buf_right;
  int *signal_A_left_ptr, *signal_A_right_ptr, *signal_B_left_ptr, *signal_B_right_ptr;

  MPI_Win win_rcv_buf_left, win_rcv_buf_right;
  MPI_Win win_signal_A_left, win_signal_A_right, win_signal_B_left, win_signal_B_right;
  int offset_left, offset_right;
  MPI_Comm comm_sm; int size_comm_sm;

/* Naming conventions                                                                */
/* Processes:                                                                        */
/*     my_rank-1                        my_rank                         my_rank+1    */
/* "left neighbor"                     "myself"                     "right neighbor" */
/*   ...    rcv_buf_right <--- snd_buf_left snd_buf_right ---> rcv_buf_left    ...   */
/*   ... snd_buf_right ---> rcv_buf_left       rcv_buf_right <--- snd_buf_left ...   */
/*                        |                                  |                       */
/*              halo-communication                 halo-communication                */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;

  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_sm);
  MPI_Comm_size(comm_sm, &size_comm_sm); 
  if (size_comm_sm != size) { printf("Not on one shared memory node \n"); MPI_Abort(MPI_COMM_WORLD, 0); }

  /* ParaStation MPI may no allow MPI_Win_allocate_shared on MPI_COMM_WORLD. Workaround: Substitute MPI_COMM_WORLD by comm_sm (on 6 lines!!!): */ 
  MPI_Win_allocate_shared((MPI_Aint)(max_length*sizeof(float)), sizeof(float), MPI_INFO_NULL, comm_sm, &rcv_buf_left,  &win_rcv_buf_left );
  MPI_Win_allocate_shared((MPI_Aint)(max_length*sizeof(float)), sizeof(float), MPI_INFO_NULL, comm_sm, &rcv_buf_right, &win_rcv_buf_right);

  MPI_Win_allocate_shared((MPI_Aint)(signal_distance*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_A_left_ptr, &win_signal_A_left);
  *signal_A_left_ptr = 0;
  MPI_Win_fence(0,win_signal_A_left);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_A_left);

  MPI_Win_allocate_shared((MPI_Aint)(signal_distance*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_A_right_ptr, &win_signal_A_right);
  *signal_A_right_ptr = 0;
  MPI_Win_fence(0,win_signal_A_right);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_A_right);

  MPI_Win_allocate_shared((MPI_Aint)(signal_distance*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_B_left_ptr, &win_signal_B_left);
  *signal_B_left_ptr = 0;
  MPI_Win_fence(0,win_signal_B_left);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_B_left);

  MPI_Win_allocate_shared((MPI_Aint)(signal_distance*sizeof(int)), sizeof(int), MPI_INFO_NULL, comm_sm,&signal_B_right_ptr, &win_signal_B_right);
  *signal_B_right_ptr = 0;
  MPI_Win_fence(0,win_signal_B_right);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_signal_B_right);

/*offset_left  is defined so that rcv_buf_left(xxx+offset_left) in process 'my_rank' is the same location as */
/*                                rcv_buf_left(xxx) in process 'left':                                       */
  offset_left  = +(left-my_rank)*max_length;

/*offset_right is defined so that rcv_buf_right(xxx+offset_right) in process 'my_rank' is the same location as */
/*                                rcv_buf_right(xxx) in process 'right':                                       */
  offset_right  = +(right-my_rank)*max_length;

  if (my_rank == 0) printf("    message size      transfertime  duplex bandwidth per process and neighbor\n");

  length = start_length;

  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_rcv_buf_left);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_rcv_buf_right);

  MPI_Win_sync(win_rcv_buf_right); /* called by the origin process */
  MPI_Win_sync(win_rcv_buf_left);  /* called by the origin process */
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_sync(win_rcv_buf_right); /* called by the target process */
  MPI_Win_sync(win_rcv_buf_left);  /* called by the target process */

  for (j = 1; j <= number_package_sizes; j++)
  { 
    
    for (i = 0; i <= number_of_messages; i++)
    {
      if(i==1) start = MPI_Wtime();

      test_value = j*1000000 + i*10000 + my_rank*10 ; mid = (length-1)/number_of_messages*i;

      snd_buf_left[0]=test_value+1  ; snd_buf_left[mid]=test_value+2  ; snd_buf_left[length-1]=test_value+3;
      snd_buf_right[0]=test_value+6 ; snd_buf_right[mid]=test_value+7 ; snd_buf_right[length-1]=test_value+8;

/*    ... The local Win_syncs are needed to sync the processor and real memory. */
/*    ... The following pair of syncs is needed that the read-write-rule is fulfilled. */
/*    rcv_buf is now exposed, i.e., the reading of previous value is finished */
      MPI_Win_sync(win_rcv_buf_right); /* called by the target process */
      MPI_Win_sync(win_rcv_buf_left); /* called by the target process */

/*    ... posting to left and right that rcv_bufs are exposed, i.e.,
          the neighbor processes are now allowed to store data into the local rcv_buf  */
      volatile_store(win_signal_A_left,  signal_A_left_ptr+signal_distance*(left-my_rank), 1);
      volatile_store(win_signal_A_right, signal_A_right_ptr+signal_distance*(right-my_rank), 1);
      while (volatile_load(win_signal_A_left,  signal_A_left_ptr)==0)  /*IDLE*/;
      while (volatile_load(win_signal_A_right, signal_A_right_ptr)==0) /*IDLE*/;
      /* while ((volatile_load(win_signal_A_left,  signal_A_left_ptr) & volatile_load(win_signal_A_right, signal_A_right_ptr))==0) -*IDLE*/;
      volatile_store(win_signal_A_left,  signal_A_left_ptr, 0);
      volatile_store(win_signal_A_right, signal_A_right_ptr, 0);
/*    ... the four syncs for the both signals can be removed, because
          the following  MPI_Win_sync(win_rcv_buf_...) will sync all sixvariables */
/*    MPI_Win_sync(win_signal_A_left); MPI_Win_sync(win_signal_A_right); */
/*    MPI_Win_sync(win_signal_B_left); MPI_Win_sync(win_signal_B_right); */

      MPI_Win_sync(win_rcv_buf_left);  /* called by the origin process */
      MPI_Win_sync(win_rcv_buf_right); /* called by the origin process */
/*    rcv_bufs are now writable by the origin process */     

/*    ...the origin process writes into the target rcv_bufs */
/*    MPI_Put(snd_buf_left,  length, MPI_FLOAT, left,  (MPI_Aint)0, length, MPI_FLOAT, win_rcv_buf_right); */
/*    MPI_Put(snd_buf_right, length, MPI_FLOAT, right, (MPI_Aint)0, length, MPI_FLOAT, win_rcv_buf_left ); */
/*      ... is substited by: */
      if (length > 1000000)
      {
         for(k=0; k<length; k++) rcv_buf_right[k+offset_left]  = snd_buf_left [k];
         for(k=0; k<length; k++) rcv_buf_left [k+offset_right] = snd_buf_right[k];
      } else {
         for(k=0; k<length; k++)
         {
            rcv_buf_right[k+offset_left]  = snd_buf_left [k];
            rcv_buf_left [k+offset_right] = snd_buf_right[k];
         }
      }

/*    ... The following pair of syncs is needed that the write-read-rule is fulfilled. */
/*    writing of rcv_buf from remote is finished */
      MPI_Win_sync(win_rcv_buf_right);  /* called by the origin process */
      MPI_Win_sync(win_rcv_buf_left);   /* called by the origin process */

/*    ... The following communication synchronizes the processors in the way */
/*        that the origin processor has finished the store */
/*        before the target processor starts to load the data.   */
/*    ... posting to left and right that rcv_bufs were stored from right and left */
      volatile_store(win_signal_B_left,  signal_B_left_ptr+signal_distance*(left-my_rank), 1);
      volatile_store(win_signal_B_right, signal_B_right_ptr+signal_distance*(right-my_rank), 1);
      while (volatile_load(win_signal_B_left,  signal_B_left_ptr)==0) /*IDLE*/;
      while (volatile_load(win_signal_B_right, signal_B_right_ptr)==0)/*IDLE*/;
      /* while ((volatile_load(win_signal_B_left,  signal_B_left_ptr) & volatile_load(win_signal_B_right, signal_B_right_ptr))==0)  -*IDLE*/;
      volatile_store(win_signal_B_left,  signal_B_left_ptr, 0);
      volatile_store(win_signal_B_right, signal_B_right_ptr, 0);
/*    ... the four syncs for the both signals can be removed, because
          the following  MPI_Win_sync(win_rcv_buf_...) will sync all sixvariables */
/*    MPI_Win_sync(win_signal_A_left); MPI_Win_sync(win_signal_A_right); */
/*    MPI_Win_sync(win_signal_B_left); MPI_Win_sync(win_signal_B_right); */

      MPI_Win_sync(win_rcv_buf_right);  /* called by the target process */
      MPI_Win_sync(win_rcv_buf_left);   /* called by the target process */
/*    rcv_buf is now locally (i.e., by the origin process) readable */ 

/*    ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process */
      test_value = j*1000000 + i*10000 + left*10  ; mid = (length-1)/number_of_messages*i;
      snd_buf_right[0]=test_value+6 ; snd_buf_right[mid]=test_value+7 ; snd_buf_right[length-1]=test_value+8;
      test_value = j*1000000 + i*10000 + right*10 ; mid = (length-1)/number_of_messages*i;
      snd_buf_left[0]=test_value+1  ; snd_buf_left[mid]=test_value+2  ; snd_buf_left[length-1]=test_value+3;
      if ((rcv_buf_left[0] != snd_buf_right[0]) || (rcv_buf_left[mid] != snd_buf_right[mid]) || 
                                                   (rcv_buf_left[length-1] != snd_buf_right[length-1])) {
         printf("%d: j=%d, i=%d --> snd_buf_right[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf_right[0], snd_buf_right[mid], snd_buf_right[length-1]);
         printf("%d:     is not identical to rcv_buf_left[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank,       mid, length-1, rcv_buf_left[0],  rcv_buf_left[mid],  rcv_buf_left[length-1]);
      }
      if ((rcv_buf_right[0] != snd_buf_left[0]) || (rcv_buf_right[mid] != snd_buf_left[mid]) ||
                                                   (rcv_buf_right[length-1] != snd_buf_left[length-1])) {
         printf("%d: j=%d, i=%d --> snd_buf_left[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf_left[0],  snd_buf_left[mid],  snd_buf_left[length-1]);
         printf("%d:     is not identical to rcv_buf_right[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank,       mid, length-1, rcv_buf_right[0], rcv_buf_right[mid], rcv_buf_right[length-1]);
      }

    }
    finish = MPI_Wtime();

    if (my_rank == 0) 
    {
      transfer_time = (finish - start) / number_of_messages;
      printf("%10i bytes %12.3f usec %13.3f MB/s\n", 
             length*(int)sizeof(float), transfer_time*1e6, 1.0e-6*2*length*sizeof(float) / transfer_time);
    }

    length = length * length_factor;
  }

  MPI_Win_unlock_all(win_rcv_buf_left);
  MPI_Win_unlock_all(win_rcv_buf_right);
  MPI_Win_unlock_all(win_signal_A_left);  MPI_Win_unlock_all(win_signal_A_right);
  MPI_Win_unlock_all(win_signal_B_left);  MPI_Win_unlock_all(win_signal_B_right);

  MPI_Win_free(&win_rcv_buf_left );
  MPI_Win_free(&win_rcv_buf_right);

  MPI_Finalize();
}
