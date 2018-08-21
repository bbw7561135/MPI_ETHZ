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
C Purpose: A program to try MPI_Issend and MPI_Recv.           C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER to_right
      PARAMETER(to_right=201)

      INTEGER ierror, my_rank, size

      INTEGER right, left

      INTEGER i, sum

      INTEGER snd_buf, rcv_buf

      INTEGER status(MPI_STATUS_SIZE)

      INTEGER request

      INTEGER(KIND=MPI_ADDRESS_KIND) iadummy


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

      sum = 0
      snd_buf = my_rank

      DO i = 1, size

         CALL MPI_ISSEND(snd_buf, 1, MPI_INTEGER, right, to_right, 
     &                   MPI_COMM_WORLD, request, ierror)

         CALL MPI_RECV(rcv_buf, 1, MPI_INTEGER, left, to_right,
     &                 MPI_COMM_WORLD, status, ierror)

         CALL MPI_WAIT(request, status, ierror)

         CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)

         snd_buf = rcv_buf
         sum = sum + rcv_buf

      END DO

      WRITE(*,*) "PE", my_rank, ": Sum =", sum

      CALL MPI_FINALIZE(ierror)

      END
