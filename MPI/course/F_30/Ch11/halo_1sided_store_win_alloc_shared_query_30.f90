PROGRAM halo

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
! Purpose: A program to meassure 1-dim halo communication      !
!          in myrank -1 and +1 directions (left and right)     !
!                                                              !
! Contents: F-Source                                           !
!                                                              !
!==============================================================!

  USE mpi_f08
  USE, INTRINSIC ::ISO_C_BINDING
  IMPLICIT NONE

  INTEGER, PARAMETER :: number_of_messages=50
  INTEGER, PARAMETER :: start_length=4
  INTEGER, PARAMETER :: length_factor=8
  INTEGER, PARAMETER :: max_length=8388608     ! ==> 2 x 32 MB per process
  INTEGER, PARAMETER :: number_package_sizes=8
! INTEGER, PARAMETER :: max_length=67108864    ! ==> 2 x 0.5 GB per process
! INTEGER, PARAMETER :: number_package_sizes=9

  INTEGER :: i, j, length, my_rank, left, right, size, test_value, mid
  INTEGER(KIND=MPI_ADDRESS_KIND) :: lb, size_of_real
  INTEGER(KIND=MPI_ADDRESS_KIND) :: iadummy
  DOUBLE PRECISION :: start, finish, transfer_time
  REAL :: snd_buf_left(max_length), snd_buf_right(max_length)
  REAL, POINTER, ASYNCHRONOUS :: rcv_buf_left(:), rcv_buf_right(:)
  REAL, POINTER, ASYNCHRONOUS :: rcv_buf_left_neighbor(:), rcv_buf_right_neighbor(:)

  TYPE(MPI_Win) :: win_rcv_buf_left, win_rcv_buf_right
  INTEGER :: disp_unit
  INTEGER(KIND=MPI_ADDRESS_KIND) :: buf_size, target_disp
  TYPE(C_PTR) :: ptr_rcv_buf_left, ptr_rcv_buf_right
  TYPE(C_PTR) :: ptr_rcv_buf_left_neighbor, ptr_rcv_buf_right_neighbor
  TYPE(MPI_Info) :: info_noncontig  
  TYPE(MPI_Comm) :: comm_sm
  INTEGER :: size_comm_sm

! Naming conventions
! Processes: 
!     my_rank-1                        my_rank                         my_rank+1
! "left neighbor"                     "myself"                     "right neighbor"
!   ...    rcv_buf_right <--- snd_buf_left snd_buf_right ---> rcv_buf_left    ...
!   ... snd_buf_right ---> rcv_buf_left       rcv_buf_right <--- snd_buf_left ...
!                        |                                  |
!              halo-communication                 halo-communication

  CALL MPI_Init()
  CALL MPI_Comm_size(MPI_COMM_WORLD, size)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  left  = mod(my_rank-1+size, size)
  right = mod(my_rank+1,      size)

  CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_sm)
  CALL MPI_Comm_size(comm_sm, size_comm_sm) 
  IF (size_comm_sm .NE. size) THEN
     write (*,*) "Not on one shared memory node"
     CALL MPI_Abort(MPI_COMM_WORLD, 0)
  END IF

  CALL MPI_Type_get_extent(MPI_REAL, lb, size_of_real) 

  CALL MPI_Info_create(info_noncontig)
  CALL MPI_Info_set(info_noncontig, "alloc_shared_noncontig", "true")
  buf_size = max_length * size_of_real
  disp_unit = size_of_real
! ...ParaStation MPI may not allow MPI_Win_allocate_shared on MPI_COMM_WORLD. Workaround: Substitute MPI_COMM_WORLD by comm_sm (on 2 lines):
  CALL MPI_Win_allocate_shared(buf_size, disp_unit, info_noncontig, MPI_COMM_WORLD, ptr_rcv_buf_left, win_rcv_buf_left)
  CALL C_F_POINTER(ptr_rcv_buf_left, rcv_buf_left, (/max_length/))
  CALL MPI_Win_allocate_shared(buf_size, disp_unit, info_noncontig, MPI_COMM_WORLD, ptr_rcv_buf_right, win_rcv_buf_right)
  CALL C_F_POINTER(ptr_rcv_buf_right, rcv_buf_right, (/max_length/))

! ... shared memory access to the rcv_buf_left, of the RIGHT neighbor process
  CALL MPI_Win_shared_query(win_rcv_buf_left, right, buf_size, disp_unit, ptr_rcv_buf_left_neighbor)
  CALL C_F_POINTER(ptr_rcv_buf_left_neighbor, rcv_buf_left_neighbor, (/max_length/))

! ... shared memory access to the rcv_buf_right, of the LEFT neighbor process
  CALL MPI_Win_shared_query(win_rcv_buf_right, left, buf_size, disp_unit, ptr_rcv_buf_right_neighbor)
  CALL C_F_POINTER(ptr_rcv_buf_right_neighbor, rcv_buf_right_neighbor, (/max_length/))

  target_disp = 0

  IF (my_rank .EQ. 0) THEN
     WRITE (*,*) "message size   transfertime    duplex bandwidth per process and neigbor"
  END IF

  length = start_length

  DO j = 1, number_package_sizes

     DO i = 0, number_of_messages

        IF (i==1) start = MPI_WTIME()

        test_value = j*1000000 + i*10000 + my_rank*10 ; mid = 1 + (length-1)/number_of_messages*i

        snd_buf_left(1)=test_value+1  ; snd_buf_left(mid)=test_value+2  ; snd_buf_left(length)=test_value+3
        snd_buf_right(1)=test_value+6 ; snd_buf_right(mid)=test_value+7 ; snd_buf_right(length)=test_value+8

!       CALL MPI_Get_address(rcv_buf_left,  iadummy)
!       CALL MPI_Get_address(rcv_buf_right, iadummy)
!       ... or with MPI-3.0 and later:
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_left)
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_right)

!       ...workaround: no assertions for shared memory - MPI-3.0 and MPI-3.1 are not clear:
        CALL MPI_Win_fence(0, win_rcv_buf_left)
        CALL MPI_Win_fence(0, win_rcv_buf_right)

!       CALL MPI_Get_address(rcv_buf_left,  iadummy)
!       CALL MPI_Get_address(rcv_buf_right, iadummy)
!       ... or with MPI-3.0 and later:
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_left)
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_right)

!       CALL MPI_Put(snd_buf_left,  length, MPI_REAL, left,  target_disp, length, MPI_REAL, win_rcv_buf_right)
!       CALL MPI_Put(snd_buf_right, length, MPI_REAL, right, target_disp, length, MPI_REAL, win_rcv_buf_left)
!       ... is substited by:
        rcv_buf_right_neighbor(1:length) = snd_buf_left(1:length)
        rcv_buf_left_neighbor(1:length)  = snd_buf_right(1:length)

!       CALL MPI_Get_address(rcv_buf_left_neighbor,  iadummy)
!       CALL MPI_Get_address(rcv_buf_right_neighbor, iadummy)
!       ... or with MPI-3.0 and later:
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_left_neighbor)
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_right_neighbor)

!       ...workaround: no assertions for shared memory - MPI-3.0 and MPI-3.1 are not clear:
        CALL MPI_Win_fence(0, win_rcv_buf_left)
        CALL MPI_Win_fence(0, win_rcv_buf_right)

!       CALL MPI_Get_address(rcv_buf_left,  iadummy)
!       CALL MPI_Get_address(rcv_buf_right, iadummy)
!       ... or with MPI-3.0 and later:
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_left)
        IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_SYNC_REG(rcv_buf_right)

!       ...snd_buf_... is used to store the values that were stored in snd_buf_... in the neighbor process
        test_value = j*1000000 + i*10000 + left*10  ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf_right(1)=test_value+6 ; snd_buf_right(mid)=test_value+7 ; snd_buf_right(length)=test_value+8
        test_value = j*1000000 + i*10000 + right*10 ; mid = 1 + (length-1)/number_of_messages*i
        snd_buf_left(1)=test_value+1  ; snd_buf_left(mid)=test_value+2  ; snd_buf_left(length)=test_value+3
        IF ((rcv_buf_left(1).NE.snd_buf_right(1)).OR.(rcv_buf_left(mid).NE.snd_buf_right(mid)).OR. &
                                                     (rcv_buf_left(length).NE.snd_buf_right(length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," -->  snd_buf_right(",1,mid,length,")=", &
                                                     snd_buf_right(1),snd_buf_right(mid),snd_buf_right(length)," in left process"
           write (*,*) my_rank,":     is not identical to rcv_buf_left(",1,mid,length,")=", &
                                                     rcv_buf_left(1),rcv_buf_left(mid),rcv_buf_left(length) 
        END IF
        IF ((rcv_buf_right(1).NE.snd_buf_left(1)).OR.(rcv_buf_right(mid).NE.snd_buf_left(mid)).OR. &
                                                     (rcv_buf_right(length).NE.snd_buf_left(length))) THEN
           write (*,*) my_rank,": j=",j," i=",i," <-- snd_buf_left(",1,mid,length,")=", &
                                                     snd_buf_left(1),snd_buf_left(mid),snd_buf_left(length)," in right process"
           write (*,*) my_rank,":     is not identical to rcv_buf_right(",1,mid,length,")=", &
                                                     rcv_buf_right(1),rcv_buf_right(mid),rcv_buf_right(length) 
        END IF

     END DO
     finish = MPI_WTIME()
     
     IF (my_rank .EQ. 0) THEN
        transfer_time = (finish - start) / (number_of_messages)
        WRITE(*,*) INT(length*size_of_real),'bytes  ', transfer_time*1e6,'usec  ', 1e-6*2*length*size_of_real/transfer_time,'MB/s'
     END IF

     length = length * length_factor
  END DO

  CALL MPI_Win_free(win_rcv_buf_left)
  CALL MPI_Win_free(win_rcv_buf_right)

  CALL MPI_Finalize()
END PROGRAM
