PROGRAM ring_derived

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
! Purpose: A program that uses derived data-types.             !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER :: ierror, my_rank, size

  INTEGER :: right, left

  INTEGER :: i

  TYPE t
     SEQUENCE
     INTEGER :: ___   ! PLEASE SUBSTITUTE ALL SKELETON CODE: ____ 
     REAL    :: ___
  END TYPE t
  TYPE(t), ASYNCHRONOUS :: snd_buf
  TYPE(t) :: rcv_buf, sum

  INTEGER(KIND=MPI_ADDRESS_KIND) :: first_var_address, second_var_address
  INTEGER :: send_recv_type

  INTEGER :: array_of_block_length(2)
  INTEGER :: array_of_types(2)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(2)

  INTEGER :: status(MPI_STATUS_SIZE)

  INTEGER :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

  right = mod(my_rank+1,      size)
  left  = mod(my_rank-1+size, size)
!     ... this SPMD-style neighbor computation with modulo has the same meaning as:
!     right = my_rank + 1
!     IF (right .EQ. size) right = 0
!     left = my_rank - 1
!     IF (left .EQ. -1) left = size-1

! Create derived datatype. 

  array_of_block_length(1) = ___
  array_of_block_length(2) = ___

  array_of_types(1) = ___
  array_of_types(2) = ___

  CALL MPI_GET_ADDRESS(snd_buf%i, first_var_address, ierror)
  CALL MPI_GET_ADDRESS(snd_buf%r, second_var_address, ierror)

  array_of_displacements(1) = 0
  array_of_displacements(2) = ___

  CALL MPI_TYPE_CREATE_STRUCT(___ ... ___, send_recv_type, ierror)
  CALL MPI_TYPE_COMMIT(___, ierror)

! ---------- original source code from MPI/course/F_20/Ch4/ring_20.f90 - PLEASE MODIFY: 
  sum = 0
  snd_buf = my_rank

  DO i = 1, size

     CALL MPI_ISSEND(snd_buf, 1, MPI_INTEGER, right, to_right, MPI_COMM_WORLD, request, ierror)

     CALL MPI_RECV(rcv_buf, 1, MPI_INTEGER, left, to_right, MPI_COMM_WORLD, status, ierror)

     CALL MPI_WAIT(request, status, ierror)

     CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)
!    ... or with MPI-3.0 and later:
!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)

     snd_buf = rcv_buf
     sum = sum + rcv_buf

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum

  CALL MPI_FINALIZE(ierror)

END PROGRAM
