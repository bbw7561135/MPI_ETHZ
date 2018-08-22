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
#include <stdlib.h>
#include <mpi.h>

#define to_right 201


int main (int argc, char *argv[])
{
  int my_world_rank, world_size, my_new_rank, new_size;
  int inner_d0, mid_d0, outer_d0, dim0, c0, ic0, mc0, oc0;
  int inner_d1, mid_d1, outer_d1, dim1, c1, ic1, mc1, oc1;
  int idim, mdim, odim, whole_size, old_rank, new_rank, *ranks;
  int dims[2], periods[2], coords[2];
  int left__coord0, left__coord1, left__rank, right_coord0, right_coord1, right_rank;
  MPI_Group world_group, new_group;
  MPI_Comm  new_comm, comm_cart;

  int snd_buf, rcv_buf;
  int right, left;
  int sum, i;
  MPI_Status  status;
  MPI_Request request;


  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

  /* Large example with 2*3 cores * 2*2 CPUs * 3*3 nodes (see slides):*/
  inner_d0=2; mid_d0=2; outer_d0=3;
  inner_d1=3; mid_d1=2; outer_d1=3;
  /* A small example with 2*1 cores * 2*2 CPUs * 3*1 nodes: */
  inner_d0=2; mid_d0=2; outer_d0=3;
  inner_d1=1; mid_d1=2; outer_d1=1;
  /* Small example with 2*2 cores * 2*1 CPUs * 1*3 nodes (see slides):*/
  inner_d0=2; mid_d0=2; outer_d0=1;
  inner_d1=2; mid_d1=1; outer_d1=3;
  dim0=inner_d0*mid_d0*outer_d0; 
  dim1=inner_d1*mid_d1*outer_d1; 
  idim=inner_d0*inner_d1; mdim=mid_d0*mid_d1; 
  odim=outer_d0*outer_d1;
  whole_size=dim0*dim1;
  ranks= malloc(whole_size*sizeof(int));
  for (oc0=0; oc0<outer_d0; oc0++) /*any sequence of the nested loop works*/
   for (mc0=0; mc0<mid_d0;   mc0++)
    for (ic0=0; ic0<inner_d0; ic0++)
     for (oc1=0; oc1<outer_d1; oc1++)
      for (mc1=0; mc1<mid_d1;   mc1++)
       for (ic1=0; ic1<inner_d1; ic1++)
       {
         old_rank =    ic1 + inner_d1*ic0
                    + (mc1 + mid_d1  *mc0)*idim
                    + (oc1 + outer_d1*oc0)*idim*mdim;
         c0 = ic0 + inner_d0*mc0 + inner_d0*mid_d0*oc0;
         c1 = ic1 + inner_d1*mc1 + inner_d1*mid_d1*oc1;
         new_rank = c1 + dim1*c0;       
         TODO: ranks[..._rank] = ..._rank;
       }

  /*only for debug-print:*/
  if (my_world_rank==0) {
    for (new_rank=0; new_rank<whole_size; new_rank++) printf("new= %3i old= %3i\n", new_rank, ranks[new_rank]);
    printf("\n");
    printf("     c1=");
    for (c1=0; c1<dim1; c1++) printf(" %3i", c1);  printf("\n");
    printf("old_rank");
    for (c1=0; c1<dim1; c1++) printf("----");  printf("\n");
    for (c0=0; c0<dim0; c0++) {
      printf("c0= %2i |", c0);
      for (c1=0; c1<dim1; c1++) printf(" %3i", ranks[c1+dim1*c0]);  printf("\n");
    }
    printf("\n");
    printf("     c1=");
    for (c1=0; c1<dim1; c1++) printf(" %3i", c1);  printf("\n");
    printf("new_rank");
    for (c1=0; c1<dim1; c1++) printf("----");  printf("\n");
    for (c0=0; c0<dim0; c0++) {
      printf("c0= %2i |", c0);
      for (c1=0; c1<dim1; c1++) printf(" %3i", c1+dim1*c0);  printf("\n");
    }
    printf("\n");
  }

  if (whole_size==world_size) { 
 
  /* Establishing new_comm with the new ranking in a array "ranks": */
    TODO: ____  
    TODO: ____  
    TODO: ____  

  /* testing the new communicator, e.g., with our ring algorithm: */
    MPI_Comm_size(new_comm, &new_size); /*should be the original one*/
    MPI_Comm_rank(new_comm, &my_new_rank);
/*-*/ /*Source code without Cartesian Virtual Toplogy */ 
/*-*/ c0 = my_new_rank / dim1; /* with 3 dims: c0 = my_new_rank / (dim1*dim2)              */
/*-*/ c1 = my_new_rank - c0*dim1;    /*        c1 = (my_new_rank - c0*(dim1*dim2) / dim2   */
/*-*/                                /*        c2 = my_new_rank - c0*(dim1*dim2) - c1*dim2 */
/*V*/ /* coordinates through a cartesian Virtual toplogy based on new_comm */
/*V*/ dims[0] = dim0; dims[1] = dim1;  periods[0]=1; periods[1]=1;
/*V*/ MPI_Cart_create(new_comm, 2, dims, periods, /*reorder=*/ 0, &comm_cart);  
/*V*/ MPI_Cart_coords(comm_cart, my_new_rank, 2, coords);
/*C*/ /* coparison of the results */
/*C*/ if (c0 != coords[0]) printf("NEWrank=%3i, ERROR in coords[0]: %2i != %2i\n", my_new_rank, c0, coords[0]);
/*C*/ if (c1 != coords[1]) printf("NEWrank=%3i, ERROR in coords[1]: %2i != %2i\n", my_new_rank, c1, coords[1]);

/* Ring in direction 0, i.e., with different c0 (and same other coord(s)) */
/* ---------------------------------------------------------------------- */

/*-*/ left__coord0 = (c0-1+dim0) % dim0;   left__rank = left__coord0*dim1+c1;
/*-*/ right_coord0 = (c0+1)      % dim0;   right_rank = right_coord0*dim1+c1;   
/*V*/ /* right=..., left=... should be substituted by one call to MPI_Cart_shift(): */
/*V*/ MPI_Cart_shift(comm_cart, /*dir=*/ 0, /*disp=*/ 1, &left, &right); 
/*C*/ if (left__rank != left ) printf("DIR=0, NEWrank=%3i, ERROR in left:  %2i != %2i\n", my_new_rank, left__rank, left);
/*C*/ if (right_rank != right) printf("DIR=0, NEWrank=%3i, ERROR in right: %2i != %2i\n", my_new_rank, right_rank, right);

    sum = 0;
    snd_buf = my_new_rank;
  
    for( i = 0; i < dim0; i++) 
    {  /* without Cartesian Virtuzal Topology, comm_cart must be substituted by new-comm */
      MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right, comm_cart, &request);
      MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right, comm_cart, &status);
      MPI_Wait(&request, &status);
      snd_buf = rcv_buf;
      sum += rcv_buf;
    }
    printf ("DIR=0, RANK world: %3i new: %3i -- coords[0]: %2i [1]: %2i -- left= %3i right= %3i -- sum= %4i\n", 
                     my_world_rank, my_new_rank,coords[0],    coords[1],   left,     right,        sum);

/* Ring in direction 1, i.e., with different c1 (and same other coord(s)) */
/* ---------------------------------------------------------------------- */

/*-*/ left__coord1 = (c1-1+dim1) % dim1;   left__rank = c0*dim1 + left__coord1;
/*-*/ right_coord1 = (c1+1)      % dim1;   right_rank = c0*dim1 + right_coord1;
/*V*/ /* right=..., left=... should be substituted by one call to MPI_Cart_shift(): */
/*V*/ MPI_Cart_shift(comm_cart, /*dir=*/ 1, /*disp=*/ 1, &left, &right); 
/*C*/ if (left__rank != left ) printf("DIR=1, NEWrank=%3i, ERROR in left:  %2i != %2i\n", my_new_rank, left__rank, left);
/*C*/ if (right_rank != right) printf("DIR=1, NEWrank=%3i, ERROR in right: %2i != %2i\n", my_new_rank, right_rank, right);

    sum = 0;
    snd_buf = my_new_rank;
  
    for( i = 0; i < dim1; i++) 
    {  /* without Cartesian Virtuzal Topology, comm_cart must be substituted by new-comm */
      MPI_Issend(&snd_buf, 1, MPI_INT, right, to_right, comm_cart, &request);
      MPI_Recv(&rcv_buf, 1, MPI_INT, left, to_right, comm_cart, &status);
      MPI_Wait(&request, &status);
      snd_buf = rcv_buf;
      sum += rcv_buf;
    }
    printf ("DIR=1, RANK world: %3i new: %3i -- coords[0]: %2i [1]: %2i -- left= %3i right= %3i -- sum= %4i\n", 
                     my_world_rank, my_new_rank,coords[0],    coords[1],   left,     right,        sum);
  }

  free(ranks); /* can be done already directly after the call to MPI_Group_incl() */

  MPI_Finalize();
}
