PROGRAM version

!==============================================================!
!                                                              !
! This file has been written as a sample solution to an        !
! exercise in a course given at the High Performance           !
! Computing Centre Stuttgart (HLRS).                           !
! It is made freely available with the understanding that      !
! every copy of this file must include this header and that    !
! HLRS take no responsibility for the use of the               !
! enclosed teaching material.                                  !
!                                                              !
! Authors: Rolf Rabenseifner (HLRS)                            !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              ! 
! Purpose: Check version of the MPI library and include file   !
!                                                              !
! Contents: Fortran-Source                                     !
!                                                              !
!==============================================================!


  USE mpi_f08
  IMPLICIT NONE

  INTEGER lib_version, lib_subversion

  CALL MPI_Init()

  CALL MPI_Get_version(lib_version, lib_subversion)
  WRITE(*,*) "Version: Library:        ", lib_version, lib_subversion
  WRITE(*,*) "Version: mpi_f08 module: ", MPI_VERSION, MPI_SUBVERSION

  CALL MPI_Finalize()

END
