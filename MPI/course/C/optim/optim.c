/****************************************************************
 *                                                              *
 * This file has been written as a sample solution to an        *
 * exercise in a course given at the High Performance           *
 * Computing Centre Stuttgart (HLRS).                           *
 * It is made freely available with the understanding that      *
 * every copy of this file must include this header and that    *
 * HLRS takes no responsibility for the use of the              *
 * enclosed teaching material.                                  *
 *                                                              *
 * Authors: Rolf Rabenseifner (HLRS)                            *
 *                                                              *
 * Contact: rabenseifner@hlrs.de                                *
 *                                                              *
 * Purpose: A program to try MPI_Issend and MPI_Recv.           *
 *                                                              *
 * Contents: C-Source                                           *
 *                                                              *
 ****************************************************************/


#include	<stdio.h>
#include	<mpi.h>

#ifndef  LOOPPAR
# define LOOPPAR       1000000
#endif

#ifndef  LOOPSERIAL
# define LOOPSERIAL     100000
#endif

#ifndef  LOOPUNBALANCE
# define LOOPUNBALANCE  100000
#endif

#ifndef  MSGSIZE
# define MSGSIZE       1000000 /*Bytes*/
#endif

#define BUFLNG MSGSIZE/sizeof(int)

double t_all, t_num, t_serial, t_comm, t_idle_num, t_idle_comm;

int my_rank, size;

/*-------------------------------------------------------------------*/

void prttime(void)
/* input arguments are the final state in the global (static) variables:  */
/* t_all, t_num, t_serial, t_comm, t_idle_num, t_idle_comm, my_rank, size */
{ double t[6], tsum[6], tmin[6], tmax[6];
  t[0]=t_all; t[1]=t_num; t[2]=t_serial; t[3]=t_comm; t[4]=t_idle_num; t[5]=t_idle_comm;
  MPI_Reduce(t, tsum, 6, MPI_DOUBLE, MPI_SUM, 0 /*=root*/, MPI_COMM_WORLD);
  MPI_Reduce(t, tmin, 6, MPI_DOUBLE, MPI_MIN, 0 /*=root*/, MPI_COMM_WORLD);
  MPI_Reduce(t, tmax, 6, MPI_DOUBLE, MPI_MAX, 0 /*=root*/, MPI_COMM_WORLD);
  if (my_rank == 0)
  { double tb; tb=tsum[0]/size/100/*percent*/;
    printf("\n");
    printf("Parallel Performance Analysis on %d MPI processes\n",size);
    printf("wall clock per process(sec)     minimum    average    maximum    max-min (over all\n");
    printf("---------------------------- ---------- ---------- ---------- ---------- processes)\n");
    printf("parallelized numerics        %10.3e %10.3e %10.3e %10.3e\n", tmin[1], tsum[1]/size, tmax[1], tmax[1]-tmin[1]);
    printf("serial numerics              %10.3e %10.3e %10.3e %10.3e\n", tmin[2], tsum[2]/size, tmax[2], tmax[2]-tmin[2]);
    printf("communication                %10.3e %10.3e %10.3e %10.3e\n", tmin[3], tsum[3]/size, tmax[3], tmax[3]-tmin[3]);
    printf("idle at end of numerics      %10.3e %10.3e %10.3e %10.3e\n", tmin[4], tsum[4]/size, tmax[4], tmax[4]-tmin[4]);
    printf("idle at end of communication %10.3e %10.3e %10.3e %10.3e\n", tmin[5], tsum[5]/size, tmax[5], tmax[5]-tmin[5]);
    printf("---------------------------- ---------- ---------- ---------- ---------- ----------\n");
    printf("total (parallel execution)   %10.3e %10.3e %10.3e\n", tmin[0], tsum[0]/size, tmax[0]);
    printf("estimated serial exec. time             %10.3e   = SerialPart+Size*ParallelPart\n",
         tsum[2]/size + tsum[1]);
    printf("estimated parallel efficience           %10.3f%%  = SerialExec/ParExec/size*100%%\n",
         (tsum[2]/size + tsum[1]) / (tsum[0]/size) / size * 100 /*percent*/);
    printf("----------------------------------------------------------------------------------\n");
    printf("\n");
    printf("wall clock per process [%%]      minimum    average    maximum    max-min (over all\n");
    printf("---------------------------- ---------- ---------- ---------- ---------- processes)\n");
    printf("parallelized numerics        %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", tmin[1]/tb, tsum[1]/size/tb, tmax[1]/tb, (tmax[1]-tmin[1])/tb);
    printf("serial numerics              %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", tmin[2]/tb, tsum[2]/size/tb, tmax[2]/tb, (tmax[2]-tmin[2])/tb);
    printf("communication                %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", tmin[3]/tb, tsum[3]/size/tb, tmax[3]/tb, (tmax[3]-tmin[3])/tb);
    printf("idle at end of numerics      %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", tmin[4]/tb, tsum[4]/size/tb, tmax[4]/tb, (tmax[4]-tmin[4])/tb);
    printf("idle at end of communication %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n", tmin[5]/tb, tsum[5]/size/tb, tmax[5]/tb, (tmax[5]-tmin[5])/tb);
    printf("---------------------------- ---------- ---------- ---------- ---------- ----------\n");
    printf("total (parallel execution)   %9.2f%% %9.2f%% %9.2f%%\n", tmin[0]/tb, tsum[0]/size/tb, tmax[0]/tb);
    printf("estimated serial exec. time             %9.2f%%  = SerialPart+Size*ParallelPart\n",
         (tsum[2]/size + tsum[1])/tb);
    printf("estimated parallel efficiency           %9.2f%%  = SerialExec/ParExec/size*100%%\n",
         (tsum[2]/size + tsum[1]) / (tsum[0]/size) / size * 100 /*percent*/);
    printf("-----------------------------------------------------------------------------------\n");
    printf("Analysis of performance loss:\n");
    printf("loss due to ...\n");
    printf("not parallelized (i.e., serial)  code   %9.2f%%  = SerialPart*(size-1)/size/ParExec\n",
         tsum[2]*(size-1)/size / tsum[0] * 100 /*percent*/);
    printf("communication                           %9.2f%%  = CommunicationPart / ParExec\n",
         tsum[3] / tsum[0] * 100 /*percent*/);
    printf("idle time at end of numerics epochs     %9.2f%%  = IdleNumericsPart  / ParExec\n",
         tsum[4] / tsum[0] * 100 /*percent*/);
    printf("idle time at end of communication epochs%9.2f%%  = IdleCommunicPart  / ParExec\n",
         tsum[5] / tsum[0] * 100 /*percent*/);
    printf("--------------------------------------- ---------- --------------------------------\n");
    printf("total loss                              %9.2f%%  = sum\n",
         (tsum[2]*(size-1)/size + tsum[3] + tsum[4] + tsum[5]) / tsum[0] * 100 /*percent*/);
    printf("approximated parallel efficiency        %9.2f%%  = 100%% - total loss\n",
         100 - (tsum[2]*(size-1)/size + tsum[3] + tsum[4] + tsum[5]) / tsum[0] * 100 /*percent*/);
    printf("-----------------------------------------------------------------------------------\n");
    printf("\n");
  }
}

/*-------------------------------------------------------------------*/

int huge (int val)
{ double s, x;
  int i;
  s = 0; x = val;

  t_serial -= MPI_Wtime();  /* PROFILING: numerics, are not parallelized, BEGIN */
  for (i=0; i<LOOPSERIAL; i++)
    { x = x / 2;
      s = s + x;
    }
  t_serial += MPI_Wtime();  /* PROFILING: numerics, are not parallelized, END   */

  /* the unbalanced part is now balanced, therefore
     - unbalanced loop is removed, and
     - this loop is now enlarged by the average number of iterations,
       see comment at unbalanced loop in nonoptim.c
  */
  for (i=0; i<(LOOPPAR + LOOPUNBALANCE/4); i++)
    { x = x / 2;
      s = s + x;
    }
  return ( (int) (s + 0.5) /*rounding*/ );
}

/*-------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
  int snd1buf[BUFLNG], rcv1buf[BUFLNG], snd2buf[BUFLNG], rcv2buf[BUFLNG];
  /* 1: send to right, receive from left */
  /* 2: send to left, receive from right */
  int right, left;
  int sum1, sum2, i;
  MPI_Status  status[2];
  MPI_Request rq[2];


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

# ifndef NO_IDLE
    MPI_Barrier(MPI_COMM_WORLD);  /* that all processes are starting measurement at nearly the same time */
# endif
  t_all=0; t_num=0; t_serial=0; t_comm=0; t_idle_num=0; t_idle_comm=0;
  t_all -= MPI_Wtime(); /* PROFILING: total execution time, BEGIN */

  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;

  sum1 = 0;  snd1buf[0] = my_rank;
  sum2 = 0;  snd2buf[0] = my_rank;

  for( i = 0; i < size; i++)
  {
#   ifndef NO_IDLE
      t_idle_num -= MPI_Wtime();  /* PROFILING: idle at end of numerics epoch, BEGIN */
      MPI_Barrier(MPI_COMM_WORLD);
      t_idle_num += MPI_Wtime();  /* PROFILING: idle at end of numerics epoch, END   */
#   endif

    t_comm -= MPI_Wtime();  /* PROFILING: communication epoch, BEGIN */
    MPI_Irecv(rcv1buf, BUFLNG, MPI_INT, left,  111, MPI_COMM_WORLD, &rq[0]);
    MPI_Send (snd1buf, BUFLNG, MPI_INT, right, 111, MPI_COMM_WORLD);

    MPI_Irecv(rcv2buf, BUFLNG, MPI_INT, right, 222, MPI_COMM_WORLD, &rq[1]);
    MPI_Send (snd2buf, BUFLNG, MPI_INT, left,  222, MPI_COMM_WORLD);
    t_comm += MPI_Wtime();  /* PROFILING: communication epoch, END   */

    MPI_Waitall(2, rq, status);

#   ifndef NO_IDLE
      t_idle_comm -= MPI_Wtime();  /* PROFILING: idle at end of communication epoch, BEGIN */
      MPI_Barrier(MPI_COMM_WORLD);
      t_idle_comm += MPI_Wtime();  /* PROFILING: idle at end of communication epoch, END   */
#   endif

    sum1 = sum1 + huge( rcv1buf[0] );
    snd1buf[0] = rcv1buf[0];
    sum2 = sum2 + huge( rcv2buf[0] );
    snd2buf[0] = rcv2buf[0];
  }

  printf ("PE%03i: Sum1 = %4i,  Sum2 = %4i\n", my_rank, sum1, sum2);

# ifndef NO_IDLE
    t_idle_num -= MPI_Wtime();  /* PROFILING: idle at end of numerics epoch, BEGIN */
    MPI_Barrier(MPI_COMM_WORLD);
    t_idle_num += MPI_Wtime();  /* PROFILING: idle at end of numerics epoch, END   */
# endif
  t_all += MPI_Wtime(); /* PROFILING: total execution time, END   */
  t_num =  t_all - (t_serial + t_comm + t_idle_num + t_idle_comm);

  prttime();

  MPI_Finalize();
}
