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
!          With Start-Post-Complete-Wait synchronization.      !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER :: my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER :: snd_buf
  INTEGER, ASYNCHRONOUS :: rcv_buf

  TYPE(MPI_Win) :: win 
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: integer_size, lb, iadummy
  INTEGER(KIND=MPI_ADDRESS_KIND) :: rcv_buf_size, target_disp

  TYPE(MPI_Group) :: grp_world, grp_left, grp_right


  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size)

  right = mod(my_rank+1,      size)
  left  = mod(my_rank-1+size, size)
!     ... this SPMD-style neighbor computation with modulo has the same meaning as:
!     right = my_rank + 1
!     IF (right .EQ. size) right = 0
!     left = my_rank - 1
!     IF (left .EQ. -1) left = size-1

! GET NEAREST NEIGHBOUR RANKS AS GROUPS

  CALL MPI_Comm_group(MPI_COMM_WORLD, grp_world)
  CALL MPI_Group_incl(grp_world, 1, (/left/),  grp_left)
  CALL MPI_Group_incl(grp_world, 1, (/right/), grp_right)
  CALL MPI_Group_free(grp_world)

! CREATE THE WINDOW.

  CALL MPI_TYPE_GET_EXTENT(MPI_INTEGER, lb, integer_size)
  rcv_buf_size = 1 * integer_size
  disp_unit = integer_size
  CALL MPI_WIN_CREATE(rcv_buf, rcv_buf_size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, win)
  target_disp = 0

  sum = 0
  snd_buf = my_rank

  DO i = 1, size

!    ... The compiler may move the read access to rcv_buf
!        in the previous loop iteration after the following 
!        1-sided MPI calls, because the compiler has no chance
!        to see, that rcv_buf will be modified by the following 
!        1-sided MPI calls.  Therefore a dummy routine must be 
!        called with rcv_buf as argument:
 
!    CALL MPI_GET_ADDRESS(rcv_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
 
!    ... Now, the compiler expects that rcv_buf was modified,
!        because the compiler cannot see that MPI_GET_ADDRESS
!        did nothing. Therefore the compiler cannot move any
!        access to rcv_buf across this dummy call. 
 
     CALL MPI_Win_post (grp_left,  MPI_MODE_NOSTORE, win)
     CALL MPI_Win_start(grp_right, 0, win)
 
     CALL MPI_Put(snd_buf, 1, MPI_INTEGER, right, target_disp, 1, MPI_INTEGER, win)
 
     CALL MPI_Win_complete(win)
     CALL MPI_Win_wait    (win)
 
!    ... The compiler has no chance to see, that rcv_buf was
!        modified. Therefore a dummy routine must be called
!        with rcv_buf as argument:
 
!    CALL MPI_GET_ADDRESS(rcv_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
 
!    ... Now, the compiler expects that rcv_buf was modified,
!        because the compiler cannot see that MPI_GET_ADDRESS
!        did nothing. Therefore the compiler will use the new
!        value on the memory, instead of some old value in a
!        register.

     snd_buf = rcv_buf
     sum = sum + rcv_buf

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum

  CALL MPI_Finalize()

END PROGRAM
