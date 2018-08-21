
/************************************************************************
 * This file has been written as a sample solution to an exercise in a 
 * course given at the Edinburgh Parallel Computing Centre. It is made
 * freely available with the understanding that every copy of this file
 * must include this header and that EPCC takes no responsibility for
 * the use of the enclosed teaching material.
 *
 * Author:      Joel Malard, Rolf Rabenseifner     
 *
 * Contact:     epcc-tec@epcc.ed.ac.uk, rabenseifner@hlrs.de
 *
 * Purpose:     A program to try out non-blocking point-to-point 
 *              communications.
 *
 * Contents:    C source code.
 *
 ************************************************************************/

#include	<stdio.h>
#include	<mpi.h>
#define tag_to_right 201

void main (int argc, char *argv[])
{
int ierror, snd_buf, my_rank, size;
int right, left;
int sum, i;

#define MAXDIMS 2 
MPI_Comm    new_comm, slice_comm;
int         dims[MAXDIMS], remain_dims[MAXDIMS],
            periods[MAXDIMS],
            reorder,
            coords[MAXDIMS];

MPI_Status  recv_status;

    MPI_Init(&argc, &argv);

    /* Get process info. */

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Set cartesian topology. */

    dims[0] = 0; 
    dims[1] = 0; 
    MPI_Dims_create(size, MAXDIMS, dims); 
    periods[0] = 1;
    periods[1] = 0;
    reorder = 1;
 
    MPI_Cart_create(MPI_COMM_WORLD,MAXDIMS,dims,periods,reorder,&new_comm);
    MPI_Comm_rank(new_comm, &my_rank);
    MPI_Cart_coords(new_comm, my_rank, MAXDIMS, coords); 
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(new_comm, remain_dims, &slice_comm); 

    /* Compute global sum. */
 
    snd_buf = my_rank;
    MPI_Allreduce(&snd_buf, &sum, 1, MPI_INT, MPI_SUM, slice_comm); 

    printf ("PE%d, x=%d, y=%d: \tSum = %d\n", 
             my_rank, coords[0], coords[1], sum);

    MPI_Finalize();

}
