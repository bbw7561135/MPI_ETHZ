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
 * Purpose: Trying MPI_Allreduce in a ring topology.             *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main (int argc, char *argv[])
{
	int my_rank, size;
	int sumA, sumB;
	MPI_Group world_group, sub_group;
	

	int my_color, key, my_new_rank;
	MPI_Comm new_comm;
	
	int ranges[1][3];
	
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	 
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	
	if (my_rank > floor((size - 1)/3.0)){
		my_color = 1;
		ranges[0][0] = floor((size - 1)/3.0) + 1;
		ranges[0][1] = size - 1;
		ranges[0][3] = 1;
	}
	else{
		my_color = 0;
		ranges[0][0] = 0;
		ranges[0][1] = floor((size - 1)/3.0);
		ranges[0][3] = 1;
	}
	// split
	key = 0;
	MPI_Comm_split(MPI_COMM_WORLD, my_color, key, & new_comm);
	MPI_Comm_rank(new_comm, &my_new_rank);
	
	/*// group
	MPI_Group_range_incl(world_group, 1, ranges, &sub_group);
	MPI_Comm_create(MPI_COMM_WORLD, sub_group, &new_comm);
	MPI_Comm_rank(new_comm, &my_new_rank);
	*/
	//sum for old ranks
	MPI_Allreduce (&my_rank, &sumA, 1, MPI_INT, MPI_SUM, new_comm);
	
	MPI_Allreduce(&my_new_rank, &sumB, 1, MPI_INT, MPI_SUM, new_comm);

	/*if (my_color ==0){
		MPI_Allreduce (&my_rank, &sumA, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	else{
		MPI_Allreduce(&my_new_rank, &sumB, 1, MPI_INT, MPI_SUM, new_comm);
	}*/
	/* Compute sum of all ranks. */
	printf ("PE world: %i, color=%i sub:%i\tSumA = %i\tSumB = %i\n", my_rank, my_color, my_new_rank, sumA, sumB);

	MPI_Finalize();
}
