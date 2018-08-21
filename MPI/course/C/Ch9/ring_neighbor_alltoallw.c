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
 * Purpose: Creating a 1-dimensional topology.                  *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

#define to_right 201
#define max_dims 1 


int main (int argc, char *argv[])
{
  int my_rank, size;
  int snd_buf, rcv_buf;
  int right, left;
  int sum, i;

  MPI_Comm    new_comm;
  int         dims[max_dims],
              periods[max_dims],
              reorder;
/*int         my_coords[max_dims]; */
  MPI_Aint    snd_displs[2], rcv_displs[2];
  int         snd_counts[2], rcv_counts[2]; 
  MPI_Datatype snd_types[2], rcv_types[2]; 

  MPI_Status  status;
  MPI_Request request;


  MPI_Init(&argc, &argv);

  /* Get process info. */
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Set cartesian topology. */
  dims[0] = size;
  periods[0] = 1;
  reorder = 1;
 
  MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, periods,
                  reorder,&new_comm);

  /* Get coords */
  MPI_Comm_rank(new_comm, &my_rank);
  /* MPI_Cart_coords(new_comm, my_rank, max_dims, my_coords); */ 

  /* Get nearest neighbour rank. */
  MPI_Cart_shift(new_comm, 0, 1, &left, &right);

  /* Compute global sum. */
 
  sum = 0;
  snd_buf = my_rank;

  rcv_counts[0] = 1;  MPI_Get_address(&rcv_buf, &rcv_displs[0]); snd_types[0] = MPI_INT;
  rcv_counts[1] = 0;  rcv_displs[1] = 0 /*unused*/;              snd_types[1] = MPI_INT;
  snd_counts[0] = 0;  snd_displs[0] = 0 /*unused*/;              rcv_types[0] = MPI_INT;
  snd_counts[1] = 1;  MPI_Get_address(&snd_buf, &snd_displs[1]); rcv_types[1] = MPI_INT;

  for( i = 0; i < size; i++) 
  {
    /* Substituted by MPI_Neighbor_alltoallw() :
    MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right,
                                new_comm, &request);
    MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right,
                              new_comm, &status);
    MPI_Wait(&request, &status);
    */    

    MPI_Neighbor_alltoallw(MPI_BOTTOM, snd_counts, snd_displs, snd_types, 
                           MPI_BOTTOM, rcv_counts, rcv_displs, rcv_types, new_comm);

    snd_buf = rcv_buf;
    sum += rcv_buf;
  }

  printf ("PE%i:\tSum = %i\n", my_rank, sum);
  /* printf ("PE%i, Coords = %i: Sum = %i\n", 
              my_rank, my_coords[0], sum); */

  MPI_Finalize();
}
