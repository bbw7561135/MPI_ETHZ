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

#define arr_lng 5
#define to_right 201
#define COUNT 2


int main (int argc, char *argv[])
{
  int my_rank, size;
  int right, left;

  struct buff{
     double f;
     int    i;
  } snd_buf[arr_lng], rcv_buf[arr_lng], sum[arr_lng];

  int i, j;

  int          array_of_blocklengths[COUNT];
  MPI_Aint     array_of_displacements[COUNT], first_var_address, second_var_address;
  MPI_Datatype array_of_types[COUNT], send_recv_type, send_recv_resized;

  MPI_Status  status;
  MPI_Request request;


  /* Get process and neighbour info. */
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
  array_of_blocklengths[0] = 1;
  array_of_blocklengths[1] = 1;

  MPI_Get_address(&snd_buf[0].f, &first_var_address);
  MPI_Get_address(&snd_buf[0].i, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = second_var_address - first_var_address;

  array_of_types[0] = MPI_DOUBLE;
  array_of_types[1] = MPI_INT;

  MPI_Type_create_struct(COUNT, array_of_blocklengths, array_of_displacements, array_of_types, &send_recv_type);
  MPI_Type_create_resized(send_recv_type, (MPI_Aint) 0, (MPI_Aint) sizeof(snd_buf[0]), &send_recv_resized);
  MPI_Type_commit(&send_recv_resized);

  /* Compute global sum. */
  for (j=0; j<arr_lng; j++)
  {
    sum[j].i = 0;            sum[j].f = 0;
    snd_buf[j].i = (j+1)*my_rank;  snd_buf[j].f = (j+1)*my_rank;  /* Step 1 = init */
    rcv_buf[j].i = -1;  rcv_buf[j].f = -1;
  }

  for( i = 0; i < size; i++) 
  { /* to check, whether the data transfer is correct, we do not transfer the last index */
    MPI_Issend(&snd_buf, arr_lng-1, send_recv_resized, right, to_right, MPI_COMM_WORLD, &request);  /* Step 2a */
    MPI_Recv(&rcv_buf, arr_lng-1, send_recv_resized, left, to_right, MPI_COMM_WORLD, &status);      /* Step 3 */
    MPI_Wait(&request, &status);                                                  /* Step 2b */
    for (j=0; j<arr_lng; j++)
    {
      snd_buf[j] = rcv_buf[j];                              /* Step 4 */
      sum[j].i += rcv_buf[j].i;  sum[j].f += rcv_buf[j].f;  /* Step 5 */
    }
  }

  if (my_rank==0)
  {
    int buf_mpi_size; 
    MPI_Aint buf_mpi_lb, buf_mpi_extent, buf_mpi_true_extent;
    MPI_Type_size(send_recv_type, &buf_mpi_size);
    MPI_Type_get_extent(send_recv_type, &buf_mpi_lb, &buf_mpi_extent);
    MPI_Type_get_true_extent(send_recv_type, &buf_mpi_lb, &buf_mpi_true_extent);
    printf("A-- MPI_Type_size:            %3d\n", buf_mpi_size);
    printf("B-- MPI_Type_get_true_extent: %3d\n", (int)buf_mpi_true_extent);
    printf("C-- MPI_Type_get_extent:      %3d\n", (int)buf_mpi_extent);
    MPI_Type_size(send_recv_resized, &buf_mpi_size);
    MPI_Type_get_extent(send_recv_resized, &buf_mpi_lb, &buf_mpi_extent);
    MPI_Type_get_true_extent(send_recv_resized, &buf_mpi_lb, &buf_mpi_true_extent);
    printf("D-- send_recv_resized:\n");
    printf("E-- MPI_Type_size:            %3d\n", buf_mpi_size);
    printf("F-- MPI_Type_get_true_extent: %3d\n", (int)buf_mpi_true_extent);
    printf("G-- MPI_Type_get_extent:      %3d\n", (int)buf_mpi_extent);
    printf("H-- sizeof:                   %3d\n", (int)sizeof(snd_buf[0]));
    MPI_Get_address(&sum[0], &first_var_address);
    MPI_Get_address(&sum[1], &second_var_address);
    printf("I-- real size is:             %3d\n", (int)(second_var_address-first_var_address));
    if (buf_mpi_extent != (second_var_address-first_var_address))
    {
      printf("J--  CAUTION:  mismatch of language type and MPI derived type: %3d != %3d\n",
             (int)buf_mpi_extent, (int)(second_var_address-first_var_address));
    }
    printf("K--\n");
    printf("L-- Expected results: for all, except the highjest j:  sum = (j+1)*(sum of all ranks)\n");
    printf("M-- For the highest j value, no data exchange is done: sum = -(number of processes)\n");
  }

  for (j=0; j<arr_lng; j++)
  {
    printf ("PE%3i j=%3i: Sum = %6i  %8.1f\n", my_rank, j, (int)sum[j].i, sum[j].f);
  }

  MPI_Finalize();
}
