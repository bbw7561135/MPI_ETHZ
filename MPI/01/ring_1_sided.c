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
	int my_rank, size;
	int snd_buf, rcv_buf;
	int right, left;
	int sum, i;

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

	sum = 0;
	snd_buf = my_rank;
	
	MPI_Win window;
	//create the window (Put)
	//MPI_Win_create(&rcv_buf, (MPI_Aint) sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	
	MPI_Win_create(&snd_buf, (MPI_Aint) sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
	
	for( i = 0; i < size; i++) 
	{
		MPI_Win_fence(MPI_MODE_NOPUT, window);
		//MPI_Put(&snd_buf, 1, MPI_INT, right, (MPI_Aint) 0, 1, MPI_INT, window);
		
		MPI_Get(&rcv_buf, 1, MPI_INT, right, (MPI_Aint) 0, 1, MPI_INT, window);
		/*MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right,
							  MPI_COMM_WORLD, &request);

		MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right,
							MPI_COMM_WORLD, &status);

		MPI_Wait(&request, &status);*/
		MPI_Win_fence(MPI_MODE_NOPUT, window);
		snd_buf = rcv_buf;
		sum += rcv_buf;
	}

	printf ("PE%i:\tSum = %i\n", my_rank, sum);

	MPI_Finalize();
}
