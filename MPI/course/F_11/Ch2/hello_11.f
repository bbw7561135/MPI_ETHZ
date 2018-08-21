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
C Purpose: A simple MPI-program printing "Hello world!"        C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"


      INTEGER ierror


      CALL MPI_INIT(ierror)

      WRITE(*,*) 'Hello world!'

      CALL MPI_FINALIZE(ierror)

      END
