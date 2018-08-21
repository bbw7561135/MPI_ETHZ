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
 * Authors: Joel Malard, Alan Simpson, (EPCC)                   *
 *          Rolf Rabenseifner          (HLRS)                   *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                *
 *                                                              *
 * Purpose: A first MPI example for this course                 *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

#include <stdio.h>
#include <mpi.h>


int main(int argc, char *argv[])
{
  int n;  double result;  // application-related data
  int my_rank, num_procs; // MPI-related data

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (my_rank == 0)   
  { // reading the application data "n" from stdin only by process 0:
    printf("Enter the number of intervals (n): \n");
    scanf("%d",&n); 
  }
  // broadcasting the content of variable "n" in process 0 
  // into variables "n" in all other processes:
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // doing some application work in each process, e.g.:
  result = 1.0 * my_rank * n;
  printf("I am process %i out of %i handling the %i's part of n=%i intervals, result=%f \n", 
                       my_rank,  num_procs,      num_procs,     n,                   result );

  if (my_rank != 0)
  { // sending some results from all processes (except 0) to process 0:
    MPI_Send(&result, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
  } else
  { // receiving all these messages and, e.g., printing them 
    int rank;
    printf("I'm proc 0: My own result is %f \n", result); 
    for (rank=1; rank<num_procs; rank++)
    {
      MPI_Recv(&result, 1, MPI_DOUBLE, rank, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("I'm proc 0: received result of process %i is %f \n", rank, result); 
    }
  }

  MPI_Finalize();
}
