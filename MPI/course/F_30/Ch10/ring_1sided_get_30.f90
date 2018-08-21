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
!          with window=snd_buf and MPI_GET to get              !
!          remote window (snd_buf) value into local rcv_buf.   !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER :: my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER, ASYNCHRONOUS :: snd_buf
  INTEGER :: rcv_buf

  TYPE(MPI_Win) :: win 
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: integer_size, lb, iadummy
  INTEGER(KIND=MPI_ADDRESS_KIND) :: rcv_buf_size, target_disp


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

! CREATE THE WINDOW.

  CALL MPI_Type_get_extent(MPI_INTEGER, lb, integer_size)
  rcv_buf_size = 1 * integer_size
  disp_unit = integer_size
  CALL MPI_Win_create(snd_buf, rcv_buf_size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, win)
  target_disp = 0

  sum = 0
  snd_buf = my_rank

  DO i = 1, size

!    ... The compiler may move the write access to snd_buf
!        in the previous loop iteration (e.g. by internally 
!        using a register instead of memory) after the following 
!        1-sided MPI calls, because the compiler has no chance
!        to see, that snd_buf will be modified by the following 
!        1-sided MPI calls.  Therefore a dummy routine must be 
!        called with rcv_buf as argument:
 
!    CALL MPI_GET_ADDRESS(snd_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
 
!    ... Now, the compiler expects that snd_buf was modified,
!        because the compiler cannot see that MPI_GET_ADDRESS
!        did nothing. Therefore the compiler cannot move any
!        access to snd_buf across this dummy call. 
 
     CALL MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE, win)
 
     CALL MPI_Get(rcv_buf, 1, MPI_INTEGER, left, target_disp, 1, MPI_INTEGER, win)
 
     CALL MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPUT + MPI_MODE_NOSUCCEED, win)
 
!    ... The compiler has no chance to see, that snd_buf was
!        read form another process between the two FENCEs. 
!        I compiler would be allowed to move the following 
!        assignment before the second FENCE. To prohibit
!        such a code move, a dummy routine must be called
!        with snd_buf as argument:
 
!    CALL MPI_GET_ADDRESS(snd_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
 
!    ... Now, the compiler expects that snd_buf was modified,
!        because the compiler cannot see that MPI_GET_ADDRESS
!        did nothing. Therefore the compiler cannot move any
!        access to snd_buf across this dummy call. 

!    CALL MPI_GET_ADDRESS(rcv_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)
!    ... This dummy call with rcv_buf in the argument list prevents
!        the following load access to rcv_buf may be moved by the
!        compiler to a position between PUT(rcv_buf...) and FENCE.

     snd_buf = rcv_buf
     sum = sum + rcv_buf

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum

  CALL MPI_Finalize()

END PROGRAM
