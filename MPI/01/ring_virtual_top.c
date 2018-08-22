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

	
	/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
	/* right = my_rank + 1;          */
	/* if (right == size) right = 0; */
	/* left = my_rank - 1;           */
	/* if (left == -1) left = size-1;*/

	sum = 0;
	//snd_buf = my_rank;
	
	
	MPI_Comm comm_cart;
	int dims;
	int coords;
	int source_rank, dest_rank;
	//int to_right = 17;
	int periods = 1;
	
	MPI_Cart_create(MPI_COMM_WORLD, 1, &size, &periods, 0, &comm_cart);
	//MPI_Cart_coords(comm_cart, my_rank, 1, &coords);
	MPI_Comm_rank(comm_cart, &my_rank);
	
	//right = (my_rank+1)      % size;
	//left  = (my_rank-1+size) % size;
	snd_buf = my_rank;
	//MPI_Cart_rank(comm_cart, &right, &source_rank);
	//MPI_Cart_rank(comm_cart, &left, &dest_rank);
	
	MPI_Cart_shift(comm_cart, 0, 1, &source_rank, &dest_rank);
	for( i = 0; i < size; i++) 
	{
		MPI_Issend(&snd_buf, 1, MPI_INT, dest_rank, to_right,
							  comm_cart, &request);

		MPI_Recv(&rcv_buf, 1, MPI_INT, source_rank, to_right,
							comm_cart, &status);

		MPI_Wait(&request, &status);

		snd_buf = rcv_buf;
		sum += rcv_buf;
	}

	printf ("PE%i:\tSum = %i\n", my_rank, sum);

	MPI_Finalize();
}
