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
! Purpose: Creating a 1-dimensional topology.                  !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER, PARAMETER :: max_dims=1

  INTEGER :: ierror, my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER, ASYNCHRONOUS :: snd_buf
  INTEGER :: rcv_buf

  INTEGER :: status(MPI_STATUS_SIZE)

  INTEGER :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy

  INTEGER :: new_comm          
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims)
! INTEGER :: coords(max_dims) 

  INTEGER(KIND=MPI_ADDRESS_KIND) :: snd_displs(2), rcv_displs(2)
  INTEGER :: snd_counts(2), rcv_counts(2)
  INTEGER :: snd_types(2), rcv_types(2)


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

! Set one-dimensional cartesian topology.

  dims(1) = size
  periods(1) = .TRUE.
  reorder = .TRUE.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm, ierror)
  CALL MPI_COMM_RANK(new_comm, my_rank, ierror)
! CALL MPI_CART_COORDS(new_comm,my_rank, max_dims, coords, ierror) 

! Get nearest neighbour ranks.

  CALL MPI_CART_SHIFT(new_comm, 0, 1, left, right, ierror)

! Compute sum.

  sum = 0
  snd_buf = my_rank

  rcv_counts(1) = 1;  CALL MPI_Get_address(rcv_buf, rcv_displs(1), ierror);  snd_types(1) = MPI_INTEGER
  rcv_counts(2) = 0;  rcv_displs(2) = 0;                                     snd_types(2) = MPI_INTEGER
  snd_counts(1) = 0;  snd_displs(1) = 0;                                     rcv_types(1) = MPI_INTEGER
  snd_counts(2) = 1;  CALL MPI_Get_address(snd_buf, snd_displs(2), ierror);  rcv_types(2) = MPI_INTEGER

  DO i = 1, size

!    ... Substituted by MPI_Neighbor_alltoallw() :
!     CALL MPI_ISSEND(snd_buf, 1, MPI_INTEGER, right, to_right, new_comm, request, ierror)
!     CALL MPI_RECV(rcv_buf, 1, MPI_INTEGER, left, to_right, new_comm, status, ierror)
!     CALL MPI_WAIT(request, status, ierror)
!     CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)
!!    ... or with MPI-3.0 and later:
!!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)

     CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)
     CALL MPI_GET_ADDRESS(rcv_buf, iadummy, ierror)
!    ... or with MPI-3.0 and later:
!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)

     CALL MPI_Neighbor_alltoallw(MPI_BOTTOM, snd_counts, snd_displs, snd_types, &
                                 MPI_BOTTOM, rcv_counts, rcv_displs, rcv_types, new_comm, ierror)

     CALL MPI_GET_ADDRESS(snd_buf, iadummy, ierror)
     CALL MPI_GET_ADDRESS(rcv_buf, iadummy, ierror)
!    ... or with MPI-3.0 and later:
!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
!    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(rcv_buf)

     snd_buf = rcv_buf
     sum = sum + rcv_buf

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum
! WRITE(*,*) "PE", my_rank, " Coord =", coords(1), ": Sum =", sum

  CALL MPI_FINALIZE(ierror)

END PROGRAM
