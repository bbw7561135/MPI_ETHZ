      MODULE global

!***************************************************************
!                                                              *
! This file has been written as a sample solution to an        *
! exercise in a course given at the High Performance           *
! Computing Centre Stuttgart (HLRS).                           *
! It is made freely available with the understanding that      *
! every copy of this file must include this header and that    *
! HLRS takes no responsibility for the use of the              *
! enclosed teaching material.                                  *
!                                                              *
! Authors: Rolf Rabenseifner (HLRS)                            *
!                                                              *
! Contact: rabenseifner@hlrs.de                                *
!                                                              *
! Purpose: A program to try MPI_Issend and MPI_Recv.           *
!                                                              *
! Contents: C-Source                                           *
!                                                              *
!***************************************************************

        IMPLICIT NONE
        DOUBLE PRECISION t_all, t_num, t_serial, t_comm
        DOUBLE PRECISION t_idle_num, t_idle_comm
        INTEGER my_rank, size

      END MODULE global

!-----------------------------------------------------------------------

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
# define MSGSIZE       1000000
#endif
!                      in Bytes

#define BUFLNG MSGSIZE/4

!-----------------------------------------------------------------------


      PROGRAM optim

        USE global

        IMPLICIT NONE

        INCLUDE 'mpif.h'

        INTEGER hugexxx

        INTEGER snd1buf(BUFLNG), rcv1buf(BUFLNG)
        INTEGER snd2buf(BUFLNG), rcv2buf(BUFLNG)
!       *** 1: send to right, receive from left
!       *** 2: send to left, receive from right
        INTEGER right, left, sum1, sum2, i
        INTEGER status(MPI_STATUS_SIZE,2), rq(2), ierror


        CALL MPI_Init(ierror, ierror)
        CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
        CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierror)

# ifndef NO_IDLE
!       *** that all processes are starting measurement at nearly the same time
        CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
# endif
        t_all=0
        t_num=0
        t_serial=0
        t_comm=0
        t_idle_num=0
        t_idle_comm=0

!       *** PROFILING: total execution time, BEGIN
        t_all = t_all - MPI_Wtime()

        right = mod(my_rank+1, size)
        left  = mod(my_rank-1+size, size)

        sum1 = 0
        snd1buf(1) = my_rank
        sum2 = 0
        snd2buf(1) = my_rank

        DO i=1,size

#         ifndef NO_IDLE
!           *** PROFILING: idle at end of numerics epoch, BEGIN
            t_idle_num = t_idle_num - MPI_Wtime()
            CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
!           *** PROFILING: idle at end of numerics epoch, END
            t_idle_num = t_idle_num + MPI_Wtime()
#         endif

!         *** PROFILING: communication epoch, BEGIN
          t_comm = t_comm - MPI_Wtime()
          CALL MPI_Irecv(rcv1buf, BUFLNG, MPI_INTEGER, left,  111,
     +                   MPI_COMM_WORLD, rq(1), ierror)
          CALL MPI_Send (snd1buf, BUFLNG, MPI_INTEGER, right, 111,
     +                   MPI_COMM_WORLD, ierror)
          CALL MPI_Irecv(rcv2buf, BUFLNG, MPI_INTEGER, right, 222,
     +                   MPI_COMM_WORLD, rq(2), ierror)
          CALL MPI_Send (snd2buf, BUFLNG, MPI_INTEGER, left,  222,
     +                   MPI_COMM_WORLD, ierror)
          CALL MPI_Waitall(2, rq, status, ierror)
!         *** PROFILING: communication epoch, END
          t_comm = t_comm + MPI_Wtime()

#         ifndef NO_IDLE
!           *** PROFILING: idle at end of communication epoch, BEGIN
            t_idle_comm = t_idle_comm - MPI_Wtime()
            CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
!           *** PROFILING: idle at end of communication epoch, END
            t_idle_comm = t_idle_comm + MPI_Wtime()
#         endif

          sum1 = sum1 + hugexxx( rcv1buf(1) )
          snd1buf(1) = rcv1buf(1)
          sum2 = sum2 + hugexxx( rcv2buf(1) )
          snd2buf(1) = rcv2buf(1)
        END DO

        print '(''PE'',I3.3,'': Sum1 = '',I4,'',  Sum2 = '',I4)',
     +                 my_rank,         sum1,             sum2

#       ifndef NO_IDLE
!         *** PROFILING: idle at end of numerics epoch, BEGIN
          t_idle_num = t_idle_num - MPI_Wtime()
          CALL MPI_Barrier(MPI_COMM_WORLD, ierror)
!         *** PROFILING: idle at end of numerics epoch, END
          t_idle_num = t_idle_num + MPI_Wtime()
#       endif
!       *** PROFILING: total execution time, END
        t_all = t_all + MPI_Wtime()

        t_num =  t_all - (t_serial + t_comm + t_idle_num + t_idle_comm)

        CALL prttime

        CALL MPI_Finalize(ierror)

      END PROGRAM

!-----------------------------------------------------------------------

      INTEGER FUNCTION hugexxx (val)

        USE global

        IMPLICIT NONE

        INCLUDE 'mpif.h'

        INTEGER val
        DOUBLE PRECISION s, x
        INTEGER i

        s = 0
        x = val

!       *** PROFILING: numerics, are not parallelized, BEGIN
        t_serial = t_serial - MPI_Wtime()

        DO i=1,LOOPSERIAL
          x = x / 2
          s = s + x
        END DO

!       *** PROFILING: numerics, are not parallelized, END
        t_serial = t_serial + MPI_Wtime()

!       *** the unbalanced part is now balanced, therefore
!       ***  - unbalanced loop is removed, and
!       ***  - this loop is now enlarged by the average number of iterations,
!       ***     see comment at unbalanced loop in nonoptim.c

        DO i=1,LOOPPAR + LOOPUNBALANCE/4
          x = x / 2
          s = s + x
        END DO

        hugexxx = s + 0.5

      END FUNCTION


!-----------------------------------------------------------------------

      SUBROUTINE prttime
!       *** input arguments are the final state in the global variables (in module global):
!       *** t_all, t_num, t_serial, t_comm, t_idle_num, t_idle_comm, my_rank, size
        USE global
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        DOUBLE PRECISION t(6), tsum(6), tmin(6), tmax(6), tb
        INTEGER ierror
        t(1)=t_all
        t(2)=t_num
        t(3)=t_serial
        t(4)=t_comm
        t(5)=t_idle_num
        t(6)=t_idle_comm

          CALL MPI_Reduce(t, tsum, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0,
     +                    MPI_COMM_WORLD, ierror)
          CALL MPI_Reduce(t, tmin, 6, MPI_DOUBLE_PRECISION, MPI_MIN, 0,
     +                    MPI_COMM_WORLD, ierror)
          CALL MPI_Reduce(t, tmax, 6, MPI_DOUBLE_PRECISION, MPI_MAX, 0,
     +                    MPI_COMM_WORLD, ierror)
        IF (my_rank .EQ. 0) THEN
!         *** tb is already divided by 100 percent:
          tb=tsum(1)/size/100
          print '('' '')'
          print '('' Parallel Performance Analysis'//
     +             ' on '',I8,'' MPI processes'')',
     +      size
          print '('' wall clock per process(sec)  '//
     +    '   minimum    average    maximum    max-min (over all'')'
          print '('' ---------------------------- '//
     +    '---------- ---------- ---------- ---------- processes)'')'
          print '('' parallelized numerics        '//
     +             ''',E10.3,1X,E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(2), tsum(2)/size, tmax(2), tmax(2)-tmin(2)
          print '('' serial numerics              '//
     +             ''',E10.3,1X,E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(3), tsum(3)/size, tmax(3), tmax(3)-tmin(3)
          print '('' communication                '//
     +             ''',E10.3,1X,E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(4), tsum(4)/size, tmax(4), tmax(4)-tmin(4)
          print '('' idle at end of numerics      '//
     +             ''',E10.3,1X,E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(5), tsum(5)/size, tmax(5), tmax(5)-tmin(5)
          print '('' idle at end of communication '//
     +             ''',E10.3,1X,E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(6), tsum(6)/size, tmax(6), tmax(6)-tmin(6)
          print '('' ---------------------------- '//
     +      '---------- ---------- ---------- ---------- ----------'')'
          print '('' total (parallel execution)   '//
     +             ''',E10.3,1X,E10.3,1X,E10.3)',
     +      tmin(1), tsum(1)/size, tmax(1)
          print '('' estimated serial exec. time  '//
     +    '           '',E10.3,''   = SerialPart+Size*ParallelPart'')',
     +      tsum(3)/size + tsum(2)
          print '('' estimated parallel efficience'//
     +    '           '',F10.3,''%  = SerialExec/ParExec/size*100%'')',
     +      (tsum(3)/size + tsum(2)) / (tsum(1)/size) / size * 100
          print '('' -----------------------------'//
     +    '-----------------------------------------------------'')'
          print '('' '')'
          print '('' wall clock per process [%]  '//
     +    '    minimum    average    maximum    max-min (over all'')'
          print '('' ---------------------------- '//
     +    '---------- ---------- ---------- ---------- processes)'')'
          print '('' parallelized numerics        '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +     tmin(2)/tb, tsum(2)/size/tb, tmax(2)/tb, (tmax(2)-tmin(2))/tb
          print '('' serial numerics              '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +     tmin(3)/tb, tsum(3)/size/tb, tmax(3)/tb, (tmax(3)-tmin(3))/tb
          print '('' communication                '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +     tmin(4)/tb, tsum(4)/size/tb, tmax(4)/tb, (tmax(4)-tmin(4))/tb
          print '('' idle at end of numerics      '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +     tmin(5)/tb, tsum(5)/size/tb, tmax(5)/tb, (tmax(5)-tmin(5))/tb
          print '('' idle at end of communication '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +     tmin(6)/tb, tsum(6)/size/tb, tmax(6)/tb, (tmax(6)-tmin(6))/tb
          print '('' ---------------------------- '//
     +    '---------- ---------- ---------- ---------- ----------'')'
          print '('' total (parallel execution)   '//
     +             ''',F9.2,''% '',F9.2,''% '',F9.2,''%'')',
     +      tmin(1)/tb, tsum(1)/size/tb, tmax(1)/tb
          print '('' estimated serial exec. time  '//
     +    '           '',F9.2,''%  = SerialPart+Size*ParallelPart'')',
     +      (tsum(3)/size + tsum(2))/tb
          print '('' estimated parallel efficiency'//
     +    '           '',F9.2,''%  = SerialExec/ParExec/size*100%'')',
     +      (tsum(3)/size + tsum(2)) / (tsum(1)/size) / size * 100
          print '('' -----------------------------'//
     +    '------------------------------------------------------'')'
          print '('' Analysis of performance loss:'')'
          print '('' loss due to ...'')'
          print '('' not parallelized (i.e., serial)  code   '//
     +             ''',F9.2,''%  = SerialPart*(size-1)/size/ParExec'')',
     +      tsum(3)*(size-1)/size / tsum(1) * 100
          print '('' communication                           '//
     +             ''',F9.2,''%  = CommunicationPart / ParExec'')',
     +      tsum(4) / tsum(1) * 100
          print '('' idle time at end of numerics epochs     '//
     +             ''',F9.2,''%  = IdleNumericsPart  / ParExec'')',
     +      tsum(5) / tsum(1) * 100
          print '('' idle time at end of communication epochs'//
     +             ''',F9.2,''%  = IdleCommunicPart  / ParExec'')',
     +      tsum(6) / tsum(1) * 100
          print '('' --------------------------------------- '//
     +             '---------- --------------------------------'')'
          print '('' total loss                              '//
     +             ''',F9.2,''%  = sum'')',
     +      (tsum(3)*(size-1)/size+tsum(4)+tsum(5)+tsum(6))/tsum(1)*100
          print '('' approximated parallel efficiency        '//
     +             ''',F9.2,''%  = 100% - total loss'')',
     +   100-(tsum(3)*(size-1)/size+tsum(4)+tsum(5)+tsum(6))/tsum(1)*100
          print '('' ----------------------------------------'//
     +             '-------------------------------------------'')'
          print '('' '')'
        END IF
      END SUBROUTINE
