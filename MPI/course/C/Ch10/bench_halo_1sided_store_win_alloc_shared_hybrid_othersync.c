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
#define number_of_repetitions 11
#define start_length 4
#define length_factor 8
#define max_length 8388608 /* ==> 2 x 32 MB per process */
#define number_package_sizes 8
/* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
/* #define number_package_sizes 9 */

int main(int argc, char *argv[])
{
  int i, j, k, r, length, my_rank, left, right, size, test_value, mid;    
  double start, finish, transfer_time[number_of_repetitions], t; 
  float snd_buf_left[max_length], snd_buf_right[max_length];
  float *rcv_buf_left, *rcv_buf_right;

  MPI_Comm subcomm_shared;
  MPI_Group subcomm_shared_group, world_group;
  int ranks_world[3], ranks_subcomm[3]; 
  int left_subcomm, myrank_subcomm, right_subcomm, size_subcomm;

  MPI_Win win_rcv_buf_left, win_rcv_buf_right;
  int offset_left, offset_right;

  MPI_Request rq[2];
  MPI_Status status_arr[2];
  int snd_dummy_left, snd_dummy_right, rcv_dummy_left, rcv_dummy_right;

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

/* #define FIXED_SHARED_SIZE 4 */
#ifdef FIXED_SHARED_SIZE
    size_subcomm = FIXED_SHARED_SIZE;
    MPI_Comm_split(MPI_COMM_WORLD, /*color:*/ my_rank/size_subcomm, /*key:*/ 0, &subcomm_shared);
#else
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, /*key:*/ 0, MPI_INFO_NULL, &subcomm_shared); 
#endif
  MPI_Comm_size(subcomm_shared, &size_subcomm);
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Comm_group(subcomm_shared, &subcomm_shared_group);
  ranks_world[0] = left; ranks_world[1] = my_rank; ranks_world[2] = right;
  MPI_Group_translate_ranks(world_group, 3, ranks_world, subcomm_shared_group, ranks_subcomm);
  left_subcomm = ranks_subcomm[0]; myrank_subcomm = ranks_subcomm[1]; right_subcomm = ranks_subcomm[2];
  /* left_subcomm and/or right_subcomm are MPI_UNDEFINED if the neighbor communication cannot use a shared memory. */
  /* printf("my_rank=%2i, size_subcomm=%2i, l/m/r=%2i/%2i/%2i, shared  l/m/r=%2i/%2i/%2i\n",
          my_rank,     size_subcomm, left, my_rank, right, left_subcomm, myrank_subcomm, right_subcomm );
     MPI_Barrier(MPI_COMM_WORLD);  */ 

  MPI_Win_allocate_shared((MPI_Aint)(max_length*sizeof(float)), sizeof(float), MPI_INFO_NULL, subcomm_shared, &rcv_buf_left,  &win_rcv_buf_left );
  MPI_Win_allocate_shared((MPI_Aint)(max_length*sizeof(float)), sizeof(float), MPI_INFO_NULL, subcomm_shared, &rcv_buf_right, &win_rcv_buf_right);
  rq[0] = MPI_REQUEST_NULL;  rq[1] = MPI_REQUEST_NULL;

/*offset_left  is defined so that rcv_buf_left(xxx+offset_left) in process 'my_rank' is the same location as */
/*                                rcv_buf_left(xxx) in process 'left':                                       */
  if (left_subcomm != MPI_UNDEFINED) offset_left  = +(left_subcomm-myrank_subcomm)*max_length;

/*offset_right is defined so that rcv_buf_right(xxx+offset_right) in process 'my_rank' is the same location as */
/*                                rcv_buf_right(xxx) in process 'right':                                       */
  if (right_subcomm != MPI_UNDEFINED) offset_right  = +(right_subcomm-myrank_subcomm)*max_length;

  if (my_rank == 0) printf("    message size      transfertime (min, 3rd min, median, 3rd max, max)                duplex bandwidth per process and neighbor (max, 3rd max, median, 3rd min, min)\n");
  if (my_rank == 0) printf("    ------------          min      3rd min       median      3rd max          max                min       3rd min        median       3rd max           max\n");

  length = start_length;

  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_rcv_buf_left);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_rcv_buf_right);

  for (j = 1; j <= number_package_sizes; j++)
  { 
    
   for (r = 0; r < number_of_repetitions; r++)
   { 
    
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i <= number_of_messages; i++)
    {
      if(i==1) start = MPI_Wtime();

      test_value = j*1000000 + i*10000 + my_rank*10 ; mid = (length-1)/number_of_messages*i;

      snd_buf_left[0]=test_value+1  ; snd_buf_left[mid]=test_value+2  ; snd_buf_left[length-1]=test_value+3;
      snd_buf_right[0]=test_value+6 ; snd_buf_right[mid]=test_value+7 ; snd_buf_right[length-1]=test_value+8;

/*      ... The local Win_syncs are needed to sync the processor and real memory on the origin-process before the processors sync */
/*          because there is no "data transfer" on the shared memory */
/*    MPI_Win_sync(win_rcv_buf_right); */
/*    MPI_Win_sync(win_rcv_buf_left); */

/*      ... tag=17: posting to left that rcv_buf_left can be stored from left / tag=23: and rcv_buf_right from right */
      if (right_subcomm != MPI_UNDEFINED) MPI_Irecv(&rcv_dummy_left,  1, MPI_INTEGER, right_subcomm, 17, subcomm_shared, &rq[0]);
      if (left_subcomm  != MPI_UNDEFINED) MPI_Irecv(&rcv_dummy_right, 1, MPI_INTEGER, left_subcomm,  23, subcomm_shared, &rq[1]);
      if (left_subcomm  != MPI_UNDEFINED) MPI_Send (&snd_dummy_right, 0, MPI_INTEGER, left_subcomm,  17, subcomm_shared);
      if (right_subcomm != MPI_UNDEFINED) MPI_Send (&snd_dummy_left,  0, MPI_INTEGER, right_subcomm, 23, subcomm_shared);
      MPI_Waitall(2, rq, status_arr);

/*      ... The local Win_syncs are needed to sync the processor and real memory on the origin-process before the processors sync */
/*          because there is no "data transfer" on the shared memory */
/*    MPI_Win_sync(win_rcv_buf_right); */
/*    MPI_Win_sync(win_rcv_buf_left); */

/*    MPI_Put(snd_buf_left,  length, MPI_FLOAT, left,  (MPI_Aint)0, length, MPI_FLOAT, win_rcv_buf_right); */
/*    MPI_Put(snd_buf_right, length, MPI_FLOAT, right, (MPI_Aint)0, length, MPI_FLOAT, win_rcv_buf_left ); */
/*      ... is substited by: */
      if (left_subcomm  != MPI_UNDEFINED) for(k=0; k<length; k++) rcv_buf_right[k+offset_left]  = snd_buf_left [k];
      if (right_subcomm != MPI_UNDEFINED) for(k=0; k<length; k++) rcv_buf_left [k+offset_right] = snd_buf_right[k];

/*      ... The local Win_syncs are needed to sync the processor and real memory on the origin-process before the processors sync */
      if (left_subcomm  != MPI_UNDEFINED) MPI_Win_sync(win_rcv_buf_right); /* because data was stored into rcv_buf_right of LEFT  neighbor */
      if (right_subcomm != MPI_UNDEFINED) MPI_Win_sync(win_rcv_buf_left);  /* because data was stored into rcv_buf_left  of RIGHT neighbor */
 
/*      ... The following communication synchronizes the processors in the way that the origin processor has finished the store */
/*           before the target processor starts to load the data.   */
/*      ... tag=17: posting to right that rcv_buf_left was stored from left / tag=23: and rcv_buf_right from right */
      if (left_subcomm  != MPI_UNDEFINED) MPI_Irecv(&rcv_dummy_left,    1, MPI_INTEGER, left_subcomm,  17, subcomm_shared, &rq[0]);
                                     else MPI_Irecv(rcv_buf_left,  length, MPI_FLOAT,   left,          17, MPI_COMM_WORLD, &rq[0]);
      if (right_subcomm != MPI_UNDEFINED) MPI_Irecv(&rcv_dummy_right,   1, MPI_INTEGER, right_subcomm, 23, subcomm_shared, &rq[1]);
                                     else MPI_Irecv(rcv_buf_right, length, MPI_FLOAT,   right,         23, MPI_COMM_WORLD, &rq[1]);
      if (right_subcomm != MPI_UNDEFINED) MPI_Send (&snd_dummy_right,   0, MPI_INTEGER, right_subcomm, 17, subcomm_shared);
                                     else MPI_Send (snd_buf_right, length, MPI_FLOAT,   right,         17, MPI_COMM_WORLD);
      if (left_subcomm  != MPI_UNDEFINED) MPI_Send (&snd_dummy_left,    0, MPI_INTEGER, left_subcomm,  23, subcomm_shared);
                                     else MPI_Send (snd_buf_left,  length, MPI_FLOAT,   left,          23, MPI_COMM_WORLD);
      MPI_Waitall(2, rq, status_arr);

/*      ... The local Win_syncs are needed to sync the processor and real memory on the target-process after the processors sync */
      if (right_subcomm != MPI_UNDEFINED) MPI_Win_sync(win_rcv_buf_right); /* because data from right should arrived in own rcv_buf_right */
      if (left_subcomm  != MPI_UNDEFINED) MPI_Win_sync(win_rcv_buf_left);  /* because data from left  should arrived in own rcv_buf_left  */

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

    } /*for i*/
    finish = MPI_Wtime();
    transfer_time[r] = (finish - start) / number_of_messages;

   } /*for r*/


   if (my_rank == 0) 
   {
      /* sorting transfer_time[...]: */
      for (r = 0; r < number_of_repetitions; r++)
        for (i = r+1; i < number_of_repetitions; i++) if (transfer_time[i]<transfer_time[r])
               { t=transfer_time[r]; transfer_time[r]=transfer_time[i]; transfer_time[i]=t; }
      /*for (r = 0; r < number_of_repetitions; r++) printf(" %12.3f", transfer_time[r]*1e6); printf("\n");*/

      printf("%10i bytes %12.3f %12.3f %12.3f %12.3f %12.3f usec %13.3f %13.3f %13.3f %13.3f %13.3f MB/s ( %5i processes, %3i processes/node)\n", 
             length*(int)sizeof(float), 
             transfer_time[1-1]*1e6,
             transfer_time[3-1]*1e6,
             transfer_time[(1+number_of_repetitions)/2-1]*1e6,
             transfer_time[number_of_repetitions-3]*1e6,
             transfer_time[number_of_repetitions-1]*1e6,
             1.0e-6*2*length*sizeof(float) / transfer_time[1-1],
             1.0e-6*2*length*sizeof(float) / transfer_time[3-1],
             1.0e-6*2*length*sizeof(float) / transfer_time[(1+number_of_repetitions)/2-1],
             1.0e-6*2*length*sizeof(float) / transfer_time[number_of_repetitions-3],
             1.0e-6*2*length*sizeof(float) / transfer_time[number_of_repetitions-1],
             size, size_subcomm);
   }

   length = length * length_factor;
  } /*for i*/

  MPI_Win_unlock_all(win_rcv_buf_left);
  MPI_Win_unlock_all(win_rcv_buf_right);

  MPI_Win_free(&win_rcv_buf_left );
  MPI_Win_free(&win_rcv_buf_right);

  MPI_Finalize();
}
