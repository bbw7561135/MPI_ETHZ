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
! Purpose: Creating a 2-dimensional topology.                  !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER, PARAMETER :: max_dims=2

  INTEGER :: my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER, ASYNCHRONOUS :: snd_buf
  INTEGER :: rcv_buf

  TYPE(MPI_Status)  :: status

  TYPE(MPI_Request) :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy

  TYPE(MPI_Comm) :: new_comm, slice_comm
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims), remain_dims(max_dims)
  INTEGER :: coords(max_dims), size_of_slice, my_rank_in_slice 


  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size)

! Set two-dimensional cartesian topology.
  dims(1) = 0       
  dims(2) = 0
  periods(1) = .TRUE.
  periods(2) = .TRUE.
  reorder = .TRUE.

  CALL MPI_Dims_create(size, max_dims, dims)
  CALL MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm)
  CALL MPI_Comm_rank(new_comm, my_rank)
  CALL MPI_Cart_coords(new_comm,my_rank, max_dims, coords) 

! Split the new-comm into slices
  remain_dims(1) = .TRUE.
  remain_dims(2) = .FALSE.
  CALL MPI_Cart_sub(new_comm, remain_dims, slice_comm)
  CALL MPI_Comm_size(slice_comm, size_of_slice)
  CALL MPI_Comm_rank(slice_comm, my_rank_in_slice)

! Get nearest neighbour ranks.

  CALL MPI_Cart_shift(slice_comm, 0, 1, left, right)

! Compute sum.

  sum = 0
  snd_buf = my_rank

  DO i = 1, size_of_slice

     CALL MPI_Issend(snd_buf, 1, MPI_INTEGER, right, to_right, slice_comm, request)

     CALL MPI_Recv(rcv_buf, 1, MPI_INTEGER, left, to_right, slice_comm, status)

     CALL MPI_Wait(request, status)

!    CALL MPI_GET_ADDRESS(snd_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)

     snd_buf = rcv_buf
     sum = sum + rcv_buf

  END DO

  WRITE(*,*) "PE",my_rank, ", coords = (",coords(1),",", coords(2), &
                 "), Slice_rank=", my_rank_in_slice, ": Sum =", sum

  CALL MPI_Finalize()

END PROGRAM
