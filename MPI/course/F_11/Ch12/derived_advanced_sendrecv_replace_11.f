      PROGRAM ring_derived

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
C Purpose: A program that uses derived data-types.             C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER to_right
      PARAMETER(to_right=201)

      INTEGER i, ierror, my_rank, size

      INTEGER right, left

      INTEGER int_buf, int_sum
      REAL real_buf, real_sum

      COMMON /send_block/ int_buf, real_buf

      INTEGER first_var_address, second_var_address
      INTEGER send_recv_type

      INTEGER array_of_block_length(2)
      INTEGER array_of_types(2)
      INTEGER array_of_displacements(2)

      INTEGER status(MPI_STATUS_SIZE)


      CALL MPI_INIT(ierror)

C Get process and neighbour info.

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

      right = mod(my_rank+1,      size)
      left  = mod(my_rank-1+size, size)
C     ... this SPMD-style neighbor computation with modulo has the same meaning as:
C     right = my_rank + 1
C     IF (right .EQ. size) right = 0
C     left = my_rank - 1
C     IF (left .EQ. -1) left = size-1

C Create derived datatype. 

      array_of_block_length(1) = 1
      array_of_block_length(2) = 1

      array_of_types(1) = MPI_INTEGER
      array_of_types(2) = MPI_REAL

      CALL MPI_ADDRESS(int_buf, first_var_address, ierror)
      CALL MPI_ADDRESS(real_buf, second_var_address, ierror)

      array_of_displacements(1) = 0
      array_of_displacements(2) = second_var_address - first_var_address

      CALL MPI_TYPE_STRUCT(2, array_of_block_length, 
     &                        array_of_displacements,
     &                        array_of_types, send_recv_type, ierror)

      CALL MPI_TYPE_COMMIT(send_recv_type, ierror)

C Initialise data
      int_sum = 0
      real_sum = 0.0

      int_buf = my_rank
      real_buf = REAL(my_rank)

C Compute sum
      DO i = 1, size

         CALL MPI_SENDRECV_REPLACE(int_buf, 1, send_recv_type, 
     &                                right, to_right, 
     &                                left, to_right, 
     &                                MPI_COMM_WORLD, status, ierror)

         int_sum = int_sum + int_buf
         real_sum = real_sum + real_buf

      END DO

      WRITE(*,*) 'PE', my_rank, ': Sum = ', int_sum, real_sum

      CALL MPI_FINALIZE(ierror)

      END
