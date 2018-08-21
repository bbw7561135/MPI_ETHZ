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
 * Authors: Rolf Rabenseifner (HLRS)                            *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                * 
 *                                                              *  
 * Purpose: A program to try MPI_Ibarrier.                      *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <mpi.h>

int main (int argc, char *argv[])
{
  int my_rank, size;
  /* in the role as sending process */
  int snd_buf_A, snd_buf_B, snd_buf_C, snd_buf_D;
  int dest, number_of_dests=0;
  int snd_finished=0;
  MPI_Request snd_rq[4];
  int total_number_of_dests; /* only for verification, should be removed in real applications */
                             /* Caution: total_number_of_dests may be less than 4, see if-statements below */
  /* in the role as receiving process */
  int rcv_buf;
  MPI_Request ib_rq;
  int ib_finished=0, rcv_flag;
  MPI_Status rcv_sts;
  int number_of_recvs=0, total_number_of_recvs; /* only for verification, should be removed in real applications */

  int round=0; /* only for verification, should be removed in real applications */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* in the role as sending process */
  dest = my_rank+1;
  if ((dest>=0) && (dest<size)) {
    snd_buf_A = 1000*my_rank + dest;  /* must not be modified until send-completion with TEST or WAIT */
    MPI_Issend(&snd_buf_A,1,MPI_INT, dest,222,MPI_COMM_WORLD, &snd_rq[number_of_dests]);
    printf("A rank: %3i - sending_: message %06i from %3i to %3i\n", my_rank, snd_buf_A, my_rank, dest);
    number_of_dests++;
  }
  dest = my_rank-2;
  if ((dest>=0) && (dest<size)) {
    snd_buf_B = 1000*my_rank + dest;  /* must not be modified until send-completion with TEST or WAIT */
    MPI_Issend(&snd_buf_B,1,MPI_INT, dest,222,MPI_COMM_WORLD, &snd_rq[number_of_dests]);
    printf("A rank: %3i - sending_: message %06i from %3i to %3i\n", my_rank, snd_buf_B, my_rank, dest);
    number_of_dests++;
  }
  dest = my_rank+4;
  if ((dest>=0) && (dest<size)) {
    snd_buf_C = 1000*my_rank + dest;  /* must not be modified until send-completion with TEST or WAIT */
    MPI_Issend(&snd_buf_C,1,MPI_INT, dest,222,MPI_COMM_WORLD, &snd_rq[number_of_dests]);
    printf("A rank: %3i - sending_: message %06i from %3i to %3i\n", my_rank, snd_buf_C, my_rank, dest);
    number_of_dests++;
  }
  dest = my_rank-7;
  if ((dest>=0) && (dest<size)) {
    snd_buf_D = 1000*my_rank + dest;  /* must not be modified until send-completion with TEST or WAIT */
    MPI_Issend(&snd_buf_D,1,MPI_INT, dest,222,MPI_COMM_WORLD, &snd_rq[number_of_dests]);
    printf("A rank: %3i - sending_: message %06i from %3i to %3i\n", my_rank, snd_buf_D, my_rank, dest);
    number_of_dests++;
  }
  while(! ib_finished) {
    /* in the role as receiving process */
    /*  MPI_IPROBE(MPI_ANY_SOURCE); If there is a message then MPI_RECV for this one message: */
    /*  Implemented in an optimized way, receiving multiple messages if available */
    rcv_flag = 1;
    /* for test reasons, the following optimizing while statement may be commented out*/ /*T*/
    while(rcv_flag) { /*it is no problem to receive as many messages as possible*/       /*T*/
      round++;                  /* only for verification, should be removed in real applications */
      if(round==10) { round=0;  /* only for verification, should be removed in real applications */
        MPI_Iprobe(MPI_ANY_SOURCE,222,MPI_COMM_WORLD, &rcv_flag, MPI_STATUS_IGNORE); 
        if(rcv_flag) {
          MPI_Recv(&rcv_buf,1,MPI_INT, MPI_ANY_SOURCE,222,MPI_COMM_WORLD, &rcv_sts);
          printf("A rank: %3i - received: message %06i from %3i to %3i\n",
                  my_rank,              rcv_buf, rcv_sts.MPI_SOURCE, my_rank);
          number_of_recvs++; /* only for verification, should be removed in real applications */
        }
      } else { rcv_flag=0; }    /* only for verification, should be removed in real applications */
    }                                                                                    /*T*/
    /* in the role as sending process */
    if(!snd_finished) { /* the following lines make only sense as long as not all MPI_ISSENDs are finished */
      /* Check whether all MPI_ISSENDs are finished */
      MPI_Testall(number_of_dests, snd_rq, &snd_finished, MPI_STATUSES_IGNORE);
      /* if all MPI_ISSENDs are finished then call MPI_IBARRIER */
      if(snd_finished) { /*i.e., the first time, i.e., only once*/
        MPI_Ibarrier(MPI_COMM_WORLD, &ib_rq);
      }
    }
    /* loop until MPI_IBARRIER finished (i.e. all processes signaled that all receives are called) */ 
    if(snd_finished) { /* the test whether the MPI_IBARRIER is finished
                          can be done only if MPI_IBARRIER is already started.
                          This ist true as soon snd_finished is true */
      MPI_Test(&ib_rq, &ib_finished, MPI_STATUS_IGNORE);
    }
  } 

  /* only for verification, should be removed in real applications: */
  MPI_Reduce(&number_of_dests, &total_number_of_dests, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&number_of_recvs, &total_number_of_recvs, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    printf("B #sends= %5i  /  #receives= %5i \n", total_number_of_dests, total_number_of_recvs);  
    if (total_number_of_dests != total_number_of_recvs) printf("C ERROR !!!! Wrong number of receives\n");
  }

  MPI_Finalize();
}
