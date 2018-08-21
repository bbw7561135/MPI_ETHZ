      PROGRAM ring

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                              C
C This file has been written as a sample solution to an        C 
C exercise in a course given at the High Performance           C
C Computing Centre Stuttgart (HLRS).                           C
C The examples are based on the examples in the MPI course of  C
C the Edinburgh Parallel Computing Centre (EPCC).              C
C It is made freely available with the understanding that      C 
C every copy of this file must include this header and that    C
C HLRS and EPCC take no responsibility for the use of the      C
C enclosed teaching material.                                  C
C                                                              C
C Authors: Joel Malard, Alan Simpson,            (EPCC)        C
C          Rolf Rabenseifner, Traugott Streicher (HLRS)        C
C                                                              C
C Contact: rabenseifner@hlrs.de                                C 
C                                                              C  
C Purpose: A program to try out one-sided communication        C
C          with window=rcv_buf and MPI_PUT to put              C
C          local snd_buf value into remote window (rcv_buf).   C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER ierror, my_rank, size

      INTEGER right, left

      INTEGER i, sum

      INTEGER snd_buf, rcv_buf

      INTEGER win 
      INTEGER disp_unit
      INTEGER (KIND=MPI_ADDRESS_KIND) integer_size, lb, iadummy
      INTEGER (KIND=MPI_ADDRESS_KIND) rcv_buf_size, target_disp


      CALL MPI_INIT(ierror)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

      right = mod(my_rank+1,      size)
      left  = mod(my_rank-1+size, size)
C     ... this SPMD-style neighbor computation with modulo has the same meaning as:
C     right = my_rank + 1
C     IF (right .EQ. size) right = 0
C     left = my_rank - 1
C     IF (left .EQ. -1) left = size-1

C CREATE THE WINDOW.

      CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, integer_size, ierror)
      rcv_buf_size = 1 * integer_size
      disp_unit = integer_size
      CALL MPI_WIN_CREATE(rcv_buf, rcv_buf_size, disp_unit,
     &                    MPI_INFO_NULL, MPI_COMM_WORLD, win, ierror)
      target_disp = 0

      sum = 0
      snd_buf = my_rank

      DO i = 1, size

C         ... The compiler may move the read access to rcv_buf
C             in the previous loop iteration after the following 
C             1-sided MPI calls, because the compiler has no chance
C             to see, that rcv_buf will be modified by the following 
C             1-sided MPI calls.  Therefore a dummy routine must be 
C             called with rcv_buf as argument:
 
          CALL MPI_GET_ADDRESS(rcv_buf, iadummy, ierror)
 
C         ... Now, the compiler expects that rcv_buf was modified,
C             because the compiler cannot see that MPI_GET_ADDRESS
C             did nothing. Therefore the compiler cannot move any
C             access to rcv_buf across this dummy call. 
 
          CALL MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE,
     &                       win, ierror)
          CALL MPI_PUT(snd_buf, 1, MPI_INTEGER, right,
     &                 target_disp, 1, MPI_INTEGER, win, ierror)
          CALL MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOPUT +
     &                       MPI_MODE_NOSUCCEED, win, ierror)
 
C         ... The compiler has no chance to see, that rcv_buf was
C             modified. Therefore a dummy routine must be called
C             with rcv_buf as argument:
 
          CALL MPI_GET_ADDRESS(rcv_buf, iadummy, ierror)
 
C         ... Now, the compiler expects that rcv_buf was modified,
C             because the compiler cannot see that MPI_GET_ADDRESS
C             did nothing. Therefore the compiler will use the new
C             value on the memory, instead of some old value in a
C             register.

          CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)
C         ... This dummy call with snd_buf in the argument list prevents
C             the following store access to snd_buf may be moved by the
C             compiler to a position between PUT(snd_buf...) and FENCE.

         snd_buf = rcv_buf
         sum = sum + rcv_buf

      END DO

      WRITE(*,*) "PE", my_rank, ": Sum =", sum

      CALL MPI_FINALIZE(ierror)

      END
