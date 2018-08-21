      PROGRAM version

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                              C
C This file has been written as a sample solution to an        C
C exercise in a course given at the High Performance           C
C Computing Centre Stuttgart (HLRS).                           C
C It is made freely available with the understanding that      C
C every copy of this file must include this header and that    C
C HLRS take no responsibility for the use of the               C
C enclosed teaching material.                                  C
C                                                              C
C Authors: Rolf Rabenseifner (HLRS)                            C
C                                                              C
C Contact: rabenseifner@hlrs.de                                C 
C                                                              C  
C Purpose: Check version of the MPI library and include file   C
C                                                              C
C Contents: Fortran-Source                                     C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE
      INCLUDE "mpif.h"

      INTEGER lib_version, lib_subversion, ierror

      CALL MPI_INIT(ierror)

      CALL MPI_GET_VERSION(lib_version, lib_subversion, ierror)
      WRITE(*,*) "Version: Library: ", lib_version, lib_subversion
      WRITE(*,*) "Version: mpif.h:  ", MPI_VERSION, MPI_SUBVERSION

      CALL MPI_FINALIZE(ierror)

      END
