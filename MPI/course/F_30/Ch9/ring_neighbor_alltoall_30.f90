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
! Purpose: Using MPI_Neighbor_alltoall for ring communication. !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201

  INTEGER, PARAMETER :: max_dims=1

  INTEGER :: my_rank, size

  INTEGER :: right, left

  INTEGER :: i, sum

  INTEGER :: snd_buf_arr(0:1)
!            snd_buf = snd_buf_arr(1)
  INTEGER :: rcv_buf_arr(0:1)
!            rcv_buf = rcv_buf_arr(0)

  TYPE(MPI_Status)  :: status

  INTEGER :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy

  TYPE(MPI_Comm) :: new_comm          
  INTEGER :: dims(max_dims)
  LOGICAL :: reorder, periods(max_dims)
! INTEGER :: coords(max_dims) 


  CALL MPI_Init()

  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size)

! Set one-dimensional cartesian topology.

  dims(1) = size
  periods(1) = .TRUE.
  reorder = .TRUE.

  CALL MPI_Cart_create(MPI_COMM_WORLD, max_dims, dims, &
                           periods, reorder, new_comm)
  CALL MPI_Comm_rank(new_comm, my_rank)
! CALL MPI_Cart_coords(new_comm,my_rank, max_dims, coords) 

! Get nearest neighbour ranks.

  CALL MPI_CART_SHIFT(new_comm, 0, 1, left, right)

! Compute sum.

  sum = 0
  snd_buf_arr(1) = my_rank

  DO i = 1, size

!      CALL MPI_Issend(snd_buf_arr(1), 1, MPI_INTEGER, right, to_right, new_comm, request)
!      CALL MPI_Recv(rcv_buf_arr(0), 1, MPI_INTEGER, left, to_right, new_comm, status)
!      CALL MPI_Wait(request, status)

     CALL MPI_Neighbor_alltoall(snd_buf_arr, 1, MPI_INTEGER, rcv_buf_arr, 1, MPI_INTEGER, new_comm)

!      CALL MPI_Get_address(snd_buf_arr, iadummy)
! !    ... or with MPI-3.0 and later:
! !    IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(snd_buf_arr)

     snd_buf_arr(1) = rcv_buf_arr(0)
     sum = sum + rcv_buf_arr(0)

  END DO

  WRITE(*,*) "PE", my_rank, ": Sum =", sum
! WRITE(*,*) "PE", my_rank, " Coord =", coords(1), ": Sum =", sum

  CALL MPI_Finalize()

END PROGRAM
