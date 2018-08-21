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
 * Purpose: A program with derived datatypes.                   *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201
#define COUNT 2


int main (int argc, char *argv[])
{
  int my_rank, size;
  int right, left;

  struct buff{
     int   ___;   /* PLEASE SUBSTITUTE ALL SKELETON CODE: ____ */
     float ___;
  } snd_buf, rcv_buf, sum;

  int i;

  int          array_of_blocklengths[COUNT];
  MPI_Aint     array_of_displacements[COUNT], first_var_address, second_var_address;
  MPI_Datatype array_of_types[COUNT], send_recv_type;

  MPI_Status  status;
  MPI_Request request;

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

  /* Set MPI datatypes for sending and receiving partial sums. */
  array_of_blocklengths[0] = ___;
  array_of_blocklengths[1] = ___;

  MPI_Get_address(&snd_buf.___, &first_var_address);
  MPI_Get_address(&snd_buf.___, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = ___;

  array_of_types[0] = ___;
  array_of_types[1] = ___;

  MPI_Type_create_struct(___ ... ___,&send_recv_type);
  MPI_Type_commit(___);

/* ---------- original source code from MPI/course/C/Ch4/ring.c - PLEASE MODIFY: */
  sum = 0;
  snd_buf = my_rank;

  for( i = 0; i < size; i++) 
  {
    MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right,
                          MPI_COMM_WORLD, &request);
   
    MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right,
                        MPI_COMM_WORLD, &status);
    
    MPI_Wait(&request, &status);
    
    snd_buf = rcv_buf;
    sum += rcv_buf;
  }

  printf ("PE%i:\tSum = %i\n", my_rank, sum);

  MPI_Finalize();
}
