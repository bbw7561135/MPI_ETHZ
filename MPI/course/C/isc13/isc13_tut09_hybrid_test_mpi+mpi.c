#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>


int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);

{
/* Test 1, isc13_tut09_hybrid_v3_all.pptx Slide 13+15 */

MPI_Aint /*IN*/ local_window_count;  double /*OUT*/  *base_ptr; 
MPI_Comm  comm_all,  comm_sm;        int  my_rank_all,  my_rank_sm,  size_sm,  disp_unit;

/*additional decl:*/
MPI_Win win_sm; int i, size_all, disp_unit_left, disp_unit_right, my_rank_world;
double *base_ptr_left, *base_ptr_right; MPI_Aint win_size_left, win_size_right; 
MPI_Comm_rank (MPI_COMM_WORLD, &my_rank_world);

/*input:*/ 
/* Test 1 and 2 is done with MPI_COMM_WORLD */
comm_all = MPI_COMM_WORLD;

local_window_count = 1000;
MPI_Comm_size (comm_all, &size_all);

MPI_Comm_rank (comm_all, &my_rank_all);
MPI_Comm_split_type (comm_all, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,  &comm_sm);
MPI_Comm_rank (comm_sm, &my_rank_sm);  MPI_Comm_size (comm_sm, &size_sm);

printf("Test 1a: [%03d %03d %02d] size_all=%3d, size_sm=%2d\n", my_rank_world, my_rank_all, my_rank_sm, size_all, size_sm); 

disp_unit = sizeof(double);  /* shared memory should contain doubles */
MPI_Win_allocate_shared (local_window_count*disp_unit, disp_unit, MPI_INFO_NULL, comm_sm,  &base_ptr, &win_sm); 

MPI_Win_fence (0, win_sm);  /*local store epoch can start*/
for (i=0; i<local_window_count; i++)  base_ptr[i] = 1000*my_rank_all + i; /* fill values into local portion */
MPI_Win_fence (0, win_sm);  /* local stores are finished, remote load epoch can start */
if (my_rank_sm > 0)         printf("Test 1b: [%03d %03d %02d] left neighbor's rightmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[-1] );
if (my_rank_sm < size_sm-1) printf("Test 1c: [%03d %03d %02d] right neighbor's leftmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[local_window_count] );

MPI_Win_free(&win_sm); 


/* Test 2, isc13_tut09_hybrid_v3-RR+GJ.pptx Slide 17+18 */

MPI_Info  info_noncontig;  
MPI_Info_create (&info_noncontig);
MPI_Info_set (info_noncontig, "alloc_shared_noncontig", "true");
MPI_Win_allocate_shared (local_window_count*disp_unit, disp_unit, info_noncontig, comm_sm, &base_ptr, &win_sm );

if (my_rank_sm > 0)         MPI_Win_shared_query (win_sm, my_rank_sm - 1, &win_size_left,  &disp_unit_left,  &base_ptr_left);
if (my_rank_sm < size_sm-1) MPI_Win_shared_query (win_sm, my_rank_sm + 1, &win_size_right, &disp_unit_right, &base_ptr_right);

printf("Test 2a: [%03d %03d %02d] disp_unit=%d, win_size_left=%lld, disp_unit_left=%d, win_size_right=%lld, disp_unit_right=%d \n", 
                                   my_rank_world, my_rank_all, my_rank_sm, disp_unit,
                                          (long long)win_size_left,  disp_unit_left, 
                                          (long long)win_size_right, disp_unit_right);
printf("Test 2b: [%03d %03d %02d] base_ptr_left=%p, base_ptr=%p, base_ptr_right=%p \n", 
                                   my_rank_world, my_rank_all, my_rank_sm, base_ptr_left, base_ptr, base_ptr_right);

MPI_Win_fence (0, win_sm);  /*local store epoch can start*/
for (i=0; i<local_window_count; i++)  base_ptr[i] = 1000*my_rank_all + i; /* fill values into local portion */

MPI_Win_fence (0, win_sm);  /* local stores are finished, remote load epoch can start */
if (my_rank_sm > 0)         printf("Test 2c: [%03d %03d %02d] left neighbor's rightmost value = %6.0lf \n", 
/*                                 my_rank_world, my_rank_all, my_rank_sm, base_ptr_left[ win_size_left/disp_unit_left - 1 ] );*/
                                   my_rank_world, my_rank_all, my_rank_sm, base_ptr_left[ local_window_count           - 1 ] );
if (my_rank_sm < size_sm-1) printf("Test 2d: [%03d %03d %02d] right neighbor's leftmost value = %6.0lf \n",
                                   my_rank_world, my_rank_all, my_rank_sm, base_ptr_right[ 0 ] );

MPI_Win_free(&win_sm); 
MPI_Comm_free(&comm_sm);



/* Test 3, isc13_tut09_hybrid_v3-RR+GJ.pptx Slide 20 */

MPI_Group group_all, group_sm; int ranges[1][3], my_rank; 

/*Input:*/
/* Test 3 is done with MPI_COMM_WORLD that should have a sequential ranking SMP nodes. */
comm_all = MPI_COMM_WORLD;
/* The size of the SMP nodes must be a multiple of size_sm: */
size_sm=2;

MPI_Comm_rank (comm_all, &my_rank);   MPI_Comm_group (comm_all, &group_all);
ranges[0][0] = (my_rank / size_sm) * size_sm;   ranges[0][1] = ranges[0][0]+size_sm-1;   ranges[0][2] = 1;
MPI_Group_range_incl (group_all, 1, ranges, &group_sm);
MPI_Comm_create (comm_all, group_sm,  &comm_sm);

MPI_Comm_rank (comm_sm, &my_rank_sm);  MPI_Comm_size (comm_sm, &size_sm);
printf("Test 3a: [%03d %02d] size_all=%3d, size_sm=%2d\n", my_rank_all, my_rank_sm, size_all, size_sm);
MPI_Group_free(&group_all);

disp_unit = sizeof(double);  /* shared memory should contain doubles */
MPI_Win_allocate_shared (local_window_count*disp_unit, disp_unit, MPI_INFO_NULL, comm_sm,  &base_ptr, &win_sm); 

MPI_Win_fence (0, win_sm);  /*local store epoch can start*/
for (i=0; i<local_window_count; i++)  base_ptr[i] = 1000*my_rank_all + i; /* fill values into local portion */
MPI_Win_fence (0, win_sm);  /* local stores are finished, remote load epoch can start */
if (my_rank_sm > 0)         printf("Test 3b: [%03d %03d %02d] left neighbor's rightmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[-1] );
if (my_rank_sm < size_sm-1) printf("Test 3c: [%03d %03d %02d] right neighbor's leftmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[local_window_count] );

MPI_Win_free(&win_sm); 
MPI_Comm_free(&comm_sm);

}



/* ------------------------------------ */

{
/* Test 4, isc13_tut09_hybrid_v3-RR+GJ.pptx Slides 132-136  */

// Input: Original communicator: MPI_Comm comm_orig; (e.g. MPI_COMM_WORLD)
//        Number of dimensions:  int      ndims = 3;
//        Global periods:        int      periods_global[] = /*e.g.*/ {1,0,1};

MPI_Comm comm_orig;
MPI_Comm comm_smp_cart, comm_nodes_cart, comm_global_cart, comm_smp_flat, comm_nodes_flat, comm_global_flat;
int i, periods_nodes[3], myrank_orig, size_smp_min, size_smp_max;
int size_smp,    myrank_smp,    mycoords_smp[3];
int size_nodes,  myrank_nodes,  mycoords_nodes[3];
int size_global, myrank_global, mycoords_global[3];

int random_key, myrank_world;
MPI_Comm_size (MPI_COMM_WORLD, &size_global);
MPI_Comm_rank (MPI_COMM_WORLD, &myrank_world);
for (i=0; i<myrank_world+1; i++) random_key = drand48()*2*size_global;
MPI_Comm_split (MPI_COMM_WORLD, 0, random_key, &comm_orig);

int periods_global[] = {1,0,1};
int ndims = 3; 
int dims_global[3];

MPI_Comm_size (comm_orig,  &size_global);
MPI_Comm_rank (comm_orig,  &myrank_orig);
// Establish a communicator on each SMP node:
MPI_Comm_split_type (comm_orig,  MPI_COMM_TYPE_SHARED,  0,  MPI_INFO_NULL,  &comm_smp_flat);
MPI_Comm_size (comm_smp_flat,  &size_smp);
int  dims_smp[] = {0,0,0};   int  periods_smp[] = {0,0,0} /*always non-period*/;
MPI_Dims_create (size_smp,  ndims,  dims_smp);
printf("Test 4-A: [%03d %03d] size_smp=%3d, ndims=%1d, dims_smp=(%3d %3d %3d)\n", 
                  myrank_world, myrank_orig, size_smp, ndims, dims_smp[0],dims_smp[1],dims_smp[2]);
MPI_Cart_create (comm_smp_flat,  ndims,  dims_smp,  periods_smp, /*reorder=*/ 1,  &comm_smp_cart);
MPI_Comm_free  (&comm_smp_flat);
MPI_Comm_rank (comm_smp_cart,  &myrank_smp);
printf("Test 4-B: [%03d %03d] myrank_smp=%3d\n", myrank_world, myrank_orig, myrank_smp);
MPI_Cart_coords (comm_smp_cart,  myrank_smp,  ndims,  mycoords_smp);
printf("Test 4-C: [%03d %03d] mycoords_smp=(%3d %3d %3d)\n", 
                  myrank_world, myrank_orig, mycoords_smp[0],mycoords_smp[1],mycoords_smp[2]);
// This source code requires that all SMP nodes have the same size. It is tested: 
MPI_Allreduce (&size_smp,  &size_smp_min, 1,  MPI_INT,  MPI_MIN,  comm_orig);
MPI_Allreduce (&size_smp,  &size_smp_max, 1,  MPI_INT,  MPI_MAX,  comm_orig);
if (size_smp_min < size_smp_max)  { printf("non-equal SMP sizes\n");  MPI_Abort (comm_orig, 1); }

// Establish the node rank. It is calculated based on the sequence of ranks in comm_orig 
// in the processes with myrank_smp == 0:
MPI_Comm_split (comm_orig, myrank_smp, 0, &comm_nodes_flat);
// Result: comm_nodes_flat combines all processes with a given myrank_smp into a separate communicator.
// Caution: The node numbering within these comm_nodes-flat may be different.
// The following source code expands the numbering from comm_nodes_flat with myrank_smp == 0
// to all node-to-node communicators:
MPI_Comm_size (comm_nodes_flat,  &size_nodes);
int dims_nodes[] =  {0,0,0};   for (i=0; i<ndims; i++) periods_nodes[i] = periods_global[i];
MPI_Dims_create (size_nodes,  ndims,  dims_nodes);
if (myrank_smp==0) {
  MPI_Cart_create (comm_nodes_flat,  ndims,  dims_nodes,  periods_nodes, 1,  &comm_nodes_cart);
  MPI_Comm_rank (comm_nodes_cart,  &myrank_nodes); 
  MPI_Comm_free  (&comm_nodes_cart); /*was needed only to calculate myrank_nodes*/
}
MPI_Comm_free (&comm_nodes_flat);
MPI_Bcast (&myrank_nodes, 1, MPI_INT, 0, comm_smp_cart);
MPI_Comm_split (comm_orig, myrank_smp, myrank_nodes, &comm_nodes_flat);
MPI_Cart_create (comm_nodes_flat,  ndims,  dims_nodes,  periods_nodes, 0,  &comm_nodes_cart);
MPI_Cart_coords (comm_nodes_cart,  myrank_nodes,  ndims,  mycoords_nodes); MPI_Comm_free (&comm_nodes_flat);

// Establish the global Cartesian communicator:
for (i=0; i<ndims; i++) { 
  dims_global[i] = dims_smp[i] * dims_nodes[i];
  mycoords_global[i] = mycoords_nodes[i] * dims_smp[i] + mycoords_smp[i];
}
myrank_global = mycoords_global[0];
for (i=1; i<ndims; i++)  { myrank_global = myrank_global * dims_global[i] + mycoords_global[i]; }
MPI_Comm_split (comm_orig,  /*color*/ 0,  myrank_global,  &comm_global_flat);
MPI_Cart_create (comm_global_flat,  ndims,  dims_global,  periods_global,  0,  &comm_global_cart); 
MPI_Comm_free (&comm_global_flat);

// Result:
//   Input was: 
//     comm_orig,  ndims,  periods_global
//   Result is:
//     comm_smp_cart,    size_smp,    myrank_smp,    dims_smp,    periods_smp,     my_coords_smp,
//     comm_nodes_cart,  size_nodes,  myrank_nodes,  dims_nodes,  periods_nodes,   my_coords_nodes,
//     comm_global_cart, size_global, myrank_global, dims_global,                  my_coords_global


int myrank_global_new;  
MPI_Comm_rank (comm_global_cart,  &myrank_global_new);

printf("Test 4a: [%03d %03d %03d] myrank_global=%3d, myrank_nodes=%3d, myrank_smp=%3d\n", 
                 myrank_world, myrank_orig, myrank_global_new, myrank_global, myrank_nodes, myrank_smp);
 
printf("Test 4b: [%03d %03d %03d] size_global=%3d, size_nodes=%3d, size_smp=%3d\n", 
                 myrank_world, myrank_orig, myrank_global_new, size_global, size_nodes, size_smp);

printf("Test 4c: [%03d %03d %03d] dims_global=(%3d %3d %3d) dims_nodes=(%3d %3d %3d) dims_smp=(%3d %3d %3d)\n", 
                 myrank_world, myrank_orig, myrank_global_new,
                                  dims_global[0],dims_global[1],dims_global[2], 
                                  dims_nodes[0],dims_nodes[1],dims_nodes[2],
                                  dims_smp[0],dims_smp[1],dims_smp[2]);

printf("Test 4d: [%03d %03d %03d] mycoords_global=(%3d %3d %3d) mycoords_nodes=(%3d %3d %3d) mycoords_smp=(%3d %3d %3d)\n", 
                 myrank_world, myrank_orig, myrank_global_new,
                                  mycoords_global[0],mycoords_global[1],mycoords_global[2], 
                                  mycoords_nodes[0],mycoords_nodes[1],mycoords_nodes[2],
                                  mycoords_smp[0],mycoords_smp[1],mycoords_smp[2]);

MPI_Comm_free(&comm_orig); /*freeing the random communicator*/
}


{
/* Test 5, isc13_tut09_hybrid_v3-RR+GJ.pptx Slides 132-136, but with ndims=1 and no Cartesian communicator  */

/*additional decl:*/
MPI_Win win_sm; int i, disp_unit, disp_unit_left, disp_unit_right;
double *base_ptr_left, *base_ptr_right; MPI_Aint win_size_left, win_size_right; 
MPI_Comm comm_orig;
MPI_Comm comm_sm, comm_nodes, comm_all;
int my_rank_orig, size_sm, my_rank_sm, size_nodes, my_rank_nodes, size_all, my_rank_all;
MPI_Aint /*IN*/ local_window_count;  double /*OUT*/  *base_ptr; 

// Input: Original communicator: MPI_Comm comm_orig; (e.g. MPI_COMM_WORLD)
local_window_count = 1000;

  int random_key, my_rank_world;
  MPI_Comm_size (MPI_COMM_WORLD, &size_all);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank_world);
  for (i=0; i<my_rank_world+1; i++) random_key = drand48()*2*size_all;
  MPI_Comm_split (MPI_COMM_WORLD, 0, random_key, &comm_orig);

  MPI_Comm_size (comm_orig,  &size_all);
  MPI_Comm_rank (comm_orig,  &my_rank_orig);

// Establish a communicator on each SMP node:
MPI_Comm_split_type (comm_orig,  MPI_COMM_TYPE_SHARED,  0,  MPI_INFO_NULL,  &comm_sm);
MPI_Comm_size (comm_sm,  &size_sm);  MPI_Comm_rank (comm_sm,  &my_rank_sm);

// Establish the node rank. It is calculated based on the sequence of ranks in comm_orig 
// in the processes with my_rank_sm == 0:
MPI_Comm_split (comm_orig, my_rank_sm, 0, &comm_nodes); 
// Result: comm_nodes combines all processes with a given my_rank_sm into a separate communicator.
//         comm_nodes is unused on processes with my_rank_sm > 0.
// Caution: The node numbering within these comm_nodes may be different.
// The following source code expands the numbering from comm_nodes with my_rank_sm == 0
// to all node-to-node communicators:
MPI_Comm_size (comm_nodes,  &size_nodes);
if (my_rank_sm==0) {
  MPI_Comm_rank (comm_nodes,  &my_rank_nodes);
  MPI_Exscan (&size_sm, &my_rank_all, 1, MPI_INT, MPI_SUM, comm_nodes); 
  // does not return value on the forst rank, therefore:
  if (my_rank_nodes == 0)  my_rank_all = 0;
}
MPI_Comm_free (&comm_nodes);
MPI_Bcast (&my_rank_nodes, 1, MPI_INT, 0, comm_sm);
MPI_Bcast (&my_rank_all, 1, MPI_INT, 0, comm_sm); my_rank_all = my_rank_all + my_rank_sm;
MPI_Comm_split (comm_orig, my_rank_sm, my_rank_nodes, &comm_nodes);
MPI_Comm_split (comm_orig,  /*color*/ 0,  my_rank_all,  &comm_all);

// Result:
//   Input was: 
//     comm_orig,  ndims,  periods_all
//   Result is:
//     comm_sm,    size_sm,    my_rank_sm,
//     comm_nodes,  size_nodes,  my_rank_nodes,
//     comm_all, size_all, my_rank_all
// The comm_nodes communicators are orthogonal to the comm_sm communicators.
// If the comm_sm communicators have different sizes with smp_size >= smp_size_min
// then only the comm_nodes on processes with my_rank_sm <= smp_size_min
// connect all the SMP nodes.


int my_rank_all_new;  
MPI_Comm_rank (comm_all,  &my_rank_all_new);

printf("Test 5a: [%03d %03d %03d] my_rank_all=%3d, my_rank_nodes=%3d, my_rank_sm=%3d\n", 
                 my_rank_world, my_rank_orig, my_rank_all_new, my_rank_all, my_rank_nodes, my_rank_sm);
 
printf("Test 5b: [%03d %03d %03d] size_all=%3d, size_nodes=%3d, size_sm=%3d\n", 
                 my_rank_world, my_rank_orig, my_rank_all_new, size_all, size_nodes, size_sm);

disp_unit = sizeof(double);  /* shared memory should contain doubles */
MPI_Win_allocate_shared (local_window_count*disp_unit, disp_unit, MPI_INFO_NULL, comm_sm,  &base_ptr, &win_sm); 

MPI_Win_fence (0, win_sm);  /*local store epoch can start*/
for (i=0; i<local_window_count; i++)  base_ptr[i] = 1000*my_rank_all + i; /* fill values into local portion */
MPI_Win_fence (0, win_sm);  /* local stores are finished, remote load epoch can start */
if (my_rank_sm > 0)         printf("Test 5c: [%03d %03d %02d] left neighbor's rightmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[-1] );
if (my_rank_sm < size_sm-1) printf("Test 5d: [%03d %03d %02d] right neighbor's leftmost value = %6.0lf \n", 
                                    my_rank_world, my_rank_all, my_rank_sm, base_ptr[local_window_count] );

MPI_Win_free(&win_sm); 

}


  MPI_Finalize();

  return 0;
}
