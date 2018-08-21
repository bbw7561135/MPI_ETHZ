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
 * Purpose: A program to try MPI_Issend and MPI_Recv.           *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201


int main (int argc, char *argv[])
{
  int my_world_rank, world_size, my_sub_rank, sub_size;
  int snd_buf, rcv_buf;
  int right, left;
  int sumA, sumB, i;
  int mycolor;
  MPI_Comm  sub_comm;

  MPI_Status  status;
  MPI_Request request;


  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

  mycolor = (my_world_rank > (world_size-1)/3);
  /* This definition of mycolor implies that the first color is 0
     --> see calculation of remote_leader below */
  MPI_Comm_split(MPI_COMM_WORLD, mycolor, 0, &sub_comm);
  MPI_Comm_size(sub_comm, &sub_size);
  MPI_Comm_rank(sub_comm, &my_sub_rank);

  right = (my_sub_rank+1)          % sub_size;
  left  = (my_sub_rank-1+sub_size) % sub_size;
/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
/* right = my_sub_rank + 1;          */
/* if (right == sub_size) right = 0; */
/* left = my_sub_rank - 1;           */
/* if (left == -1) left = sub_size-1;*/

  sumA = 0;
  snd_buf = my_world_rank;
  for( i = 0; i < sub_size; i++) 
  {
    MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right,
                          sub_comm, &request);
    MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right,
                        sub_comm, &status);
    MPI_Wait(&request, &status);
    snd_buf = rcv_buf;
    sumA += rcv_buf;
  }

  sumB = 0;
  snd_buf = my_sub_rank;
  for( i = 0; i < sub_size; i++) 
  {
    MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right,
                          sub_comm, &request);
    MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right,
                        sub_comm, &status);
    MPI_Wait(&request, &status);
    snd_buf = rcv_buf;
    sumB += rcv_buf;
  }

  printf ("PE world:%3i, color=%i sub:%3i, SumA=%3i, SumB=%3i\n", 
          my_world_rank, mycolor, my_sub_rank, sumA, sumB);

  MPI_Finalize();
}
