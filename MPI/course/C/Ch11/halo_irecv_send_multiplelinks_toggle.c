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
 * Purpose: A program to meassure 1-dim halo communication      *
 *          in myrank -1 and +1 directions (left and right)     *
 * Ammendment on June 14-15, 2018 by Rolf Rabenseifner (HLRS)   *
 * for the slides on hybrid programming MPIX_2018-HY-S_v02.pptx *
 *          Added the different schemes with several small rings*
 *          within spilt communicators:                         *
 *           - rings within a CPU                               *
 *           - rings from one CPU to other CPUs within a node   *
 *           - rings accross the nodes, using only the forst CPU*
 *           - rings accross the nodes, using all CPUs          *
 * Toggle:  snd_buf and rcv_buf are exchanged after each        *
 *          communication step to prohibit wrong reporting      *
 *          of the CPU-to-_CPU communication.                   *
 *          This benchmark still reports sending out from       *
 *          warm caches (i.e., filled caches).                  *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/

/* Commands:
 *
 * - to run it on the Cray XC40 system hazelhen.hww.de at HLRS in Stuttgart:
 *    qsub -I -lnodes=4:ppn=24,walltime=0:24:00 -q test
 *    cc -o halo_irecv_send_multiplelinks.exe halo_irecv_send_multiplelinks.c
 *    time aprun -n 96 ./halo_irecv_send_multiplelinks.exe > ./halo_irecv_send_multiplelinks_prot4.txt
 *
 * - to convert the human readable protocol into a Excel readable protocol used with the provided xlsx files:
 *     sed -e 's/\([0-9]\)\.\([0-9]\)/\1,\2/g' -e 's/ \( *[0-9,][0-9,]* *\) /;\1;/g' halo_irecv_send_multiplelinks_prot4.txt > halo_irecv_send_multiplelinks_prot4_comma.txt
 *             ^__ this part is to convert the     ^__ this part is to add ";"-separators
 *                 english numbers nn.nnn into         around the numbers
 *                 german numbers  nn,nnn
 *     vi halo_irecv_send_multiplelinks_prot4_comma.txt  --> :set fileformat=dos  --> :wq
 *
 * - read in the data into the xlsx table on A1 (e.g., in halo_irecv_send_multiplelinks_prot4.xlsx):
 *     the data from, e.g., halo_irecv_send_multiplelinks_prot4_comma.txt
 *     must be read in as ";"-separated data
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define number_of_messages 50
// #define start_length 4
// #define length_factor 8
// #define max_length 8388608 /* ==> 2 x 32 MB per process */
// #define number_package_sizes 8
// /* #define max_length 67108864    */ /* ==> 2 x 0.5 GB per process */
// /* #define number_package_sizes 9 */

#define start_length 2
#define length_factor 2
#define max_length 8388608 /* ==> 2 x 32 MB per process */
#define number_package_sizes 23 

#define num_schemes 4
#define num_cores_per_cpu 12
#define num_cpus_per_node 2

int main(int argc, char *argv[])
{
  int i, j, length, my_rank, left, right, size, test_value, mid;    
  double start, finish, transfer_time, max_transfer_time; 
  MPI_Request rq[2];
  MPI_Status status_arr[2];
  /* float snd_buf_left[max_length], snd_buf_right[max_length]; */
  /* float rcv_buf_left[max_length], rcv_buf_right[max_length]; */
  float *snd_buf_left, *snd_buf_right;
  float *rcv_buf_left, *rcv_buf_right;
  snd_buf_left  = malloc(max_length*sizeof(float));
  snd_buf_right = malloc(max_length*sizeof(float));
  rcv_buf_left  = malloc(max_length*sizeof(float));
  rcv_buf_right = malloc(max_length*sizeof(float));

  int world_my_rank, world_size, scheme_num, active_cores_per_cpu, active;
  int core_num, cpu_num, node_num, color, active_cores_per_node;
  MPI_Comm comm, active_comm;

/* Naming conventions                                                                */
/* Processes:                                                                        */
/*     my_rank-1                        my_rank                         my_rank+1    */
/* "left neighbor"                     "myself"                     "right neighbor" */
/*   ...    rcv_buf_right <--- snd_buf_left snd_buf_right ---> rcv_buf_left    ...   */
/*   ... snd_buf_right ---> rcv_buf_left       rcv_buf_right <--- snd_buf_left ...   */
/*                        |                                  |                       */
/*              halo-communication                 halo-communication                */

  MPI_Init(&argc, &argv);
 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  core_num = world_my_rank % num_cores_per_cpu;
  cpu_num  = (world_my_rank / num_cores_per_cpu) % num_cpus_per_node;
  node_num = world_my_rank / num_cores_per_cpu / num_cpus_per_node; 

  if (world_my_rank == 0) {
    printf("This benchmark runs on %2i nodes, %1i CPUs per node and %2i cores per CPU, i.e., in total on %3i cores (=MPI processes).\n", 
           world_size/num_cpus_per_node/num_cores_per_cpu, num_cpus_per_node, num_cores_per_cpu, world_size);
  }
 
  for (scheme_num=0; scheme_num < num_schemes; scheme_num++)
  {
 
    if (world_my_rank == 0) {
      printf("\nNEW scheme: %1i ", scheme_num); 
      if (scheme_num==0)
      { printf("intra-CPU:\n   duplex ring communication between all cores of a CPU\n   (executed in parallel in all CPUs)\n");
      } else if (scheme_num==1)
      { printf("intra-node:\n   duplex ring communication between cores with same core number within a node\n   (executed in parallel for all cores and all nodes)\n");
      } else if (scheme_num==2)
      { printf("inter-node (only by 1st CPU):\n   duplex ring communication between cores with same core+CPU number between the nodes\n   (executed in parallel only in by the 1st CPU in each node)\n");
      } else if (scheme_num==3)
      { printf("inter-node (by all CPU per node):\n   duplex ring communication between cores with same core+CPU number between the nodes\n   (executed in parallel by all CPU in each node)\n");
      } else
      { printf("Dummy scheme, do not use!\n");
      }
    }
   
    for (active_cores_per_cpu = (scheme_num==0 ? 2 : 1); active_cores_per_cpu <= num_cores_per_cpu ; active_cores_per_cpu++)
    {
      if (scheme_num==0)
      {
         if (core_num < active_cores_per_cpu)
         { color = world_my_rank / num_cores_per_cpu /*= number of cpu in world */;
           active_cores_per_node = active_cores_per_cpu * num_cpus_per_node;
         }else
         { color = MPI_UNDEFINED;
         }      
         MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);
   
      } else if (scheme_num==1)
      {
         if (core_num < active_cores_per_cpu)
         { color = node_num * num_cores_per_cpu + core_num /* number of core in first CPU per node */;
           active_cores_per_node = active_cores_per_cpu * num_cpus_per_node;
         }else
         { color = MPI_UNDEFINED;
         }      
         MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);
   
      } else if (scheme_num==2)
      {
         if ((core_num < active_cores_per_cpu) && (cpu_num == 0))
         { color = core_num /*= number of core in 1st CPU */;
           active_cores_per_node = active_cores_per_cpu /* only one CPU is active !!!!!!! */;
         }else
         { color = MPI_UNDEFINED;
         }      
         MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);
   
      } else if (scheme_num==3)
      {
         if (core_num < active_cores_per_cpu)
         { color = world_my_rank % (num_cores_per_cpu *num_cpus_per_node) /*= number of core in node */;
           active_cores_per_node = active_cores_per_cpu * num_cpus_per_node;
         }else
         { color = MPI_UNDEFINED;
         }      
         MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);
   
      } else 
      {
         MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &comm);
         if (world_my_rank == 0) {
            printf("This scheme %1i should not be used. It is only a dummy.\n", scheme_num); 
         }
      }
 
      if (world_my_rank == 0) {
        printf("scheme: %1i ; %2i number of active cores per CPU ; %2i number of active cores per node \n", 
                  scheme_num,  active_cores_per_cpu,                                                   active_cores_per_node); 
      }
  
       active = comm != MPI_COMM_NULL;
       MPI_Comm_split(MPI_COMM_WORLD, active, 0, &active_comm);
       if (comm != MPI_COMM_NULL)
       {

  /* ====== Begin of execution of a ring on comm ======= */
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);
  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;

  if (world_my_rank == 0) printf("    message size      transfertime  duplex bandwidth per process and neighbor (comm_size= %2i)\n", size);

  length = start_length;

  for (j = 1; j <= number_package_sizes; j++)
  { 
    
    for (i = 0; i <= number_of_messages; i++)
    {
      if(i==1) start = MPI_Wtime();

      test_value = j*1000000 + i*10000 + my_rank*10 ; mid = (length-1)/number_of_messages*i;

      snd_buf_left[0]=test_value+1  ; snd_buf_left[mid]=test_value+2  ; snd_buf_left[length-1]=test_value+3;
      snd_buf_right[0]=test_value+6 ; snd_buf_right[mid]=test_value+7 ; snd_buf_right[length-1]=test_value+8;

      MPI_Irecv(rcv_buf_right, length, MPI_FLOAT, right, 17, comm, &rq[0]);
      MPI_Irecv(rcv_buf_left,  length, MPI_FLOAT, left,  23, comm, &rq[1]);

      MPI_Send(snd_buf_left,  length, MPI_FLOAT, left,  17, comm);
      MPI_Send(snd_buf_right, length, MPI_FLOAT, right, 23, comm);
    
      MPI_Waitall(2, rq, status_arr);

/*    ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process */
      test_value = j*1000000 + i*10000 + left*10  ; mid = (length-1)/number_of_messages*i;
      snd_buf_right[0]=test_value+6 ; snd_buf_right[mid]=test_value+7 ; snd_buf_right[length-1]=test_value+8;
      test_value = j*1000000 + i*10000 + right*10 ; mid = (length-1)/number_of_messages*i;
      snd_buf_left[0]=test_value+1  ; snd_buf_left[mid]=test_value+2  ; snd_buf_left[length-1]=test_value+3;
      if ((rcv_buf_left[0] != snd_buf_right[0]) || (rcv_buf_left[mid] != snd_buf_right[mid]) || 
                                                   (rcv_buf_left[length-1] != snd_buf_right[length-1])) {
         printf("%d: j=%d, i=%d --> snd_buf_right[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf_right[0], snd_buf_right[mid], snd_buf_right[length-1]);
         printf("%d:     is not identical to rcv_buf_left[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank,       mid, length-1, rcv_buf_left[0],  rcv_buf_left[mid],  rcv_buf_left[length-1]);
      }
      if ((rcv_buf_right[0] != snd_buf_left[0]) || (rcv_buf_right[mid] != snd_buf_left[mid]) ||
                                                   (rcv_buf_right[length-1] != snd_buf_left[length-1])) {
         printf("%d: j=%d, i=%d --> snd_buf_left[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank, j, i, mid, length-1, snd_buf_left[0],  snd_buf_left[mid],  snd_buf_left[length-1]);
         printf("%d:     is not identical to rcv_buf_right[0,%d,%d]=(%f,%f,%f)\n",
                    my_rank,       mid, length-1, rcv_buf_right[0], rcv_buf_right[mid], rcv_buf_right[length-1]);
      }

      { /*toggle:*/ float *tmp;
        tmp=snd_buf_left;  snd_buf_left  = rcv_buf_right; rcv_buf_right = tmp;
        tmp=snd_buf_right; snd_buf_right = rcv_buf_left;  rcv_buf_left  = tmp;
      } 

    }
    finish = MPI_Wtime();

    transfer_time = (finish - start) / number_of_messages;
    MPI_Allreduce(&transfer_time, &max_transfer_time, 1, MPI_DOUBLE, MPI_MAX, active_comm);
    if (world_my_rank == 0) 
    {
     if(active)
     {
      printf("%10i bytes rank-0: %12.3f usec %13.3f MB/s  max: %12.3f usec  min: %13.3f MB/s  accumulated-per-node: %13.3f MB/s \n", 
             length*(int)sizeof(float), transfer_time*1e6,     1.0e-6*2*length*sizeof(float) / transfer_time,
                                        max_transfer_time*1e6, 1.0e-6*2*length*sizeof(float) / max_transfer_time,
                                                               1.0e-6*2*length*sizeof(float) / max_transfer_time * active_cores_per_node);
     } else
     {
      printf("UNEXPECTED: world_my_rank == 0 is inactive \n");
     }
    }

    length = length * length_factor;
  }
  /* ====== End of execution of a ring on comm ======= */

       } /* endif comm != MPI_COMM_NULL */
    } /* for active_cores_per_cpu */
  } /* end for scheme_num */

  MPI_Finalize();
}
