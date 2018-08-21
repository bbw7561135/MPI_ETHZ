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

  USE mpi_f08

  IMPLICIT NONE

  INTEGER, PARAMETER :: arr_lng=5

  INTEGER, PARAMETER :: to_right=201

  INTEGER :: my_rank, size

  INTEGER :: right, left

  INTEGER :: i, j

! TYPE, BIND(C) :: t
  TYPE t
     SEQUENCE
     DOUBLE PRECISION :: r
     INTEGER :: i
  END TYPE t
  TYPE(t), ASYNCHRONOUS :: snd_buf(arr_lng)
  TYPE(t) :: rcv_buf(arr_lng), sum(arr_lng)

  INTEGER(KIND=MPI_ADDRESS_KIND) :: first_var_address, second_var_address
  TYPE(MPI_Datatype) :: send_recv_type

  INTEGER :: array_of_block_length(2)
  TYPE(MPI_Datatype) :: array_of_types(2)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(2)

  TYPE(MPI_Status) :: status

  TYPE(MPI_Request) :: request

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy

  INTEGER :: buf_mpi_size; 
  INTEGER(KIND=MPI_ADDRESS_KIND) :: buf_mpi_lb, buf_mpi_extent, buf_mpi_true_extent;


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

! Create derived datatype. 

  array_of_block_length(1) = 1
  array_of_block_length(2) = 1

  array_of_types(1) = MPI_DOUBLE_PRECISION
  array_of_types(2) = MPI_INTEGER

  CALL MPI_Get_address(snd_buf(1)%r, first_var_address)
  CALL MPI_Get_address(snd_buf(1)%i, second_var_address)

  array_of_displacements(1) = 0
  array_of_displacements(2) = second_var_address - first_var_address

  CALL MPI_Type_create_struct(2, array_of_block_length, array_of_displacements, array_of_types, send_recv_type)
  CALL MPI_Type_commit(send_recv_type)

  do j=1,arr_lng
    sum(j)%i = 0;              sum(j)%r = 0
    snd_buf(j)%i = j*my_rank;  snd_buf(j)%r = j*REAL(my_rank)
    rcv_buf(j)%i = -1;         rcv_buf(j)%r = -1
  end do 

  DO i = 1, size

!    ... to check, whether the data transfer is correct, we do not transfer the last index
     CALL MPI_Issend(snd_buf, arr_lng-1, send_recv_type, right, to_right, MPI_COMM_WORLD, request)

     CALL MPI_Recv(rcv_buf, arr_lng-1, send_recv_type, left, to_right, MPI_COMM_WORLD, status)

     CALL MPI_Wait(request, status)

!    CALL MPI_GET_ADDRESS(snd_buf, iadummy)
!    ... should be substituted as soon as possible by:
     IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)

     snd_buf = rcv_buf
     do j=1,arr_lng
       sum(j)%i = sum(j)%i + rcv_buf(j)%i
       sum(j)%r = sum(j)%r + rcv_buf(j)%r
     end do 

  END DO

  if (my_rank==0) THEN
    CALL MPI_Type_size(send_recv_type, buf_mpi_size)
    CALL MPI_Type_get_extent(send_recv_type, buf_mpi_lb, buf_mpi_extent)
    CALL MPI_Type_get_true_extent(send_recv_type, buf_mpi_lb, buf_mpi_true_extent)
    WRITE(*,"(a,i3)") "A-- MPI_Type_size:            ", INT(buf_mpi_size)
    WRITE(*,"(a,i3)") "B-- MPI_Type_get_true_extent: ", INT(buf_mpi_true_extent)
    WRITE(*,"(a,i3)") "C-- MPI_Type_get_extent:      ", INT(buf_mpi_extent)
    CALL MPI_GET_ADDRESS(sum(1), first_var_address)
    CALL MPI_GET_ADDRESS(sum(2), second_var_address)
    WRITE(*,"(a,i3)") "D-- real size is:             ", INT(second_var_address-first_var_address)
    if (buf_mpi_extent .NE. second_var_address-first_var_address) then
      WRITE(*,"(a,i3,a,i3)") "E-- CAUTION:  mismatch of language type and MPI derived type: ", &
                        INT(buf_mpi_extent), " != ", INT(second_var_address-first_var_address)
    end if
    WRITE(*,"(a)") "F--"
    WRITE(*,"(a)") "G-- Expected results: for all, except the highjest j:  sum = (j+1)*(sum of all ranks)"
    WRITE(*,"(a)") "H-- For the highest j value, no data exchange is done: sum = -(number of processes)"
  ENDIF

  do j=1,arr_lng
    WRITE(*,"(a,i3,a,i3,a,i8,a,f10.1)") "PE", my_rank, " j=",j,":  Sum(j)%i =", sum(j)%i, "  Sum(j)%r =", real(sum(j)%r)
  end do 

  CALL MPI_Finalize()

END PROGRAM
