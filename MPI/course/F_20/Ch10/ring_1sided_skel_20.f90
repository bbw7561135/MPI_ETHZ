PROGRAM ring

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
! Purpose: A program to try out one-sided communication        !
!          with window=rcv_buf and MPI_PUT to put              !
!          local snd_buf value into remote window (rcv_buf).   !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER :: ierror, my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER ___ :: snd_buf  ! PLEASE SUBSTITUTE ALL SKELETON CODE: ____ 
  INTEGER ___ :: rcv_buf  ! Here, you should decide, which buffer should be ASYNCHRONOUS

  INTEGER :: win 
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: integer_size, lb, iadummy
  INTEGER(KIND=MPI_ADDRESS_KIND) :: buf_size, target_disp


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

! CREATE THE WINDOW.

  target_disp = 0  ! This "long" integer zero is needed in the call to MPI_PUT  
  CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, integer_size, ierror)
  buf_size = 1 * integer_size
  disp_unit = integer_size
  CALL MPI_WIN_CREATE(___ ... ___, win, ierror)

! ---------- original source code from MPI/course/F_20/Ch4/ring_20.f90 - PLEASE MODIFY: 
  sum = 0
  snd_buf = my_rank

  DO i = 1, size

!    The following skeleton block of MPI calls must be substituted by one-sided RMA,
!    surrounded by RMA synchronization and appropriate MPI_F_SYNC_REG (or alternatives):

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
