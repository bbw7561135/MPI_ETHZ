        PROGRAM hello

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
C Purpose: A program to try MPI_Comm_size and MPI_Comm_rank.   C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"


      INTEGER ierror, rank, size


      CALL MPI_INIT(ierror)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

      IF (rank .EQ. 0) THEN
         WRITE(*,*) 'Hello world!'
      END IF  

      WRITE(*,*) 'I am process', rank, ' out of', size  

      CALL MPI_FINALIZE(ierror)
   
      END
