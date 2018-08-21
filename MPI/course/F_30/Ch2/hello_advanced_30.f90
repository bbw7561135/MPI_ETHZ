PROGRAM hello

!==============================================================!
!                                                              !
! This file has been written as a sample solution to an        !
! exercise in a course given at the High Performance           !
! Computing Centre Stuttgart (HLRS).                           !
! The examples are based on the examples in the MPI course of  !
! the Edinburgh Parallel Computing Centre (EPCC).              !
! It is made freely available with the understanding that      !
! every copy of this file must include this header and that    !
! HLRS and EPCC take no responsibility for the use of the      !
! enclosed teaching material.                                  !
!                                                              !
! Authors: Joel Malard, Alan Simpson,            (EPCC)        !
!          Rolf Rabenseifner, Traugott Streicher (HLRS)        !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              !
! Purpose: A program to try MPI_Comm_size and MPI_Comm_rank.   !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE


  INTEGER rank, size


  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank)

  CALL MPI_Comm_size(MPI_COMM_WORLD, size)

  IF (rank .EQ. 0) THEN
     WRITE(*,*) 'Hello world!'
  END IF  

  WRITE(*,*) 'I am process', rank, ' out of', size  

  CALL MPI_Finalize()
   
END PROGRAM
