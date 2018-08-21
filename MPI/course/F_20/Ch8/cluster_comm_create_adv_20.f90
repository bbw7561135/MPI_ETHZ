PROGRAM renumber

!**************************************************************!
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
! Purpose: A program to try MPI_Issend and MPI_Recv.           !
!                                                              !
! Contents: C-Source                                           !
!                                                              !
!**************************************************************!


  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: to_right=201


  INTEGER :: my_world_rank, world_size, my_new_rank, new_size, ierror
  INTEGER :: inner_d0, mid_d0, outer_d0, dim0, c0, ic0, mc0, oc0
  INTEGER :: inner_d1, mid_d1, outer_d1, dim1, c1, ic1, mc1, oc1
  INTEGER :: idim, mdim, odim, whole_size, old_rank, new_rank
  INTEGER, ALLOCATABLE :: ranks(:)
  INTEGER :: dims(0:1), coords(0:1)
  LOGICAL :: periods(0:1)
  INTEGER :: left__coord0, left__coord1, left__rank, right_coord0, right_coord1, right_rank
  INTEGER :: world_group, new_group
  INTEGER :: new_comm, comm_cart

  INTEGER :: snd_buf, rcv_buf
  INTEGER :: right, left
  INTEGER :: sum, i
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: request


  CALL MPI_INIT(ierror)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_world_rank, ierror)

  ! Large example with 2*3 cores * 2*2 CPUs * 3*3 nodes (see slides):
  inner_d0=2; mid_d0=2; outer_d0=3
  inner_d1=3; mid_d1=2; outer_d1=3
  ! A small example with 2*1 cores * 2*2 CPUs * 3*1 nodes:
  inner_d0=2; mid_d0=2; outer_d0=3
  inner_d1=1; mid_d1=2; outer_d1=1
  ! Small example with 2*2 cores * 2*1 CPUs * 1*3 nodes (see slides):
  inner_d0=2; mid_d0=2; outer_d0=1
  inner_d1=2; mid_d1=1; outer_d1=3

  dim0=inner_d0*mid_d0*outer_d0
  dim1=inner_d1*mid_d1*outer_d1
  idim=inner_d0*inner_d1
  mdim=mid_d0*mid_d1
  odim=outer_d0*outer_d1
  whole_size=dim0*dim1
  ALLOCATE(ranks(0:whole_size-1))
  DO oc0=0, outer_d0-1  ! any sequence of the nested loop works
   DO mc0=0, mid_d0-1   
    DO ic0=0, inner_d0-1 
     DO oc1=0, outer_d1-1 
      DO mc1=0, mid_d1-1   
       DO ic1=0, inner_d1-1 
         old_rank =    ic1 + inner_d1*ic0             &
          &         + (mc1 + mid_d1  *mc0)*idim       &
          &         + (oc1 + outer_d1*oc0)*idim*mdim
         c0 = ic0 + inner_d0*mc0 + inner_d0*mid_d0*oc0
         c1 = ic1 + inner_d1*mc1 + inner_d1*mid_d1*oc1
         new_rank = c1 + dim1*c0
         ranks(new_rank) = old_rank
       END DO
      END DO
     END DO
    END DO
   END DO
  END DO

  !only for debug-print:
  IF (my_world_rank==0) THEN
    DO new_rank=0, whole_size-1
      WRITE (*,'(''new= '',I3,'' old= '',I3)') new_rank, ranks(new_rank)
    END DO
    WRITE(*,*)
    WRITE(*,'(''     c1='',1000('' '',I3))') (c1, c1=0,dim1-1)
    WRITE(*,'(''old_rank'',1000A)') ('----', c1=0,dim1-1)
    DO c0=0, dim0-1
       WRITE(*,'(''c0= '',I2,'' |'',1000('' '',I3))') c0, ( ranks(c1+dim1*c0), c1=0,dim1-1)
    END DO
    WRITE(*,*)
    WRITE(*,'(''     c1='',1000('' '',I3))') (c1, c1=0,dim1-1)
    WRITE(*,'(''new_rank'',1000A)') ('----', c1=0,dim1-1)
    DO c0=0, dim0-1
       WRITE(*,'(''c0= '',I2,'' |'',1000('' '',I3))') c0, ( c1+dim1*c0, c1=0,dim1-1)
    END DO
    WRITE(*,*)
  END IF

  IF (whole_size==world_size) THEN
 
  ! Establishing new_comm with the new ranking in a array "ranks":
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror)
    CALL MPI_GROUP_INCL(world_group, world_size, ranks, new_group, ierror)
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD, new_group, new_comm, ierror)

  ! testing the new communicator, e.g., with our ring algorithm:
    CALL MPI_COMM_SIZE(new_comm, new_size, ierror)  ! should be the original one
    CALL MPI_COMM_RANK(new_comm, my_new_rank, ierror)
    ! Source code without Cartesian Virtual Toplogy
    c0 = my_new_rank / dim1  ! with 3 dims: c0 = my_new_rank / (dim1*dim2)
    c1 = my_new_rank - c0*dim1     !        c1 = (my_new_rank - c0*(dim1*dim2) / dim2
                                   !        c2 = my_new_rank - c0*(dim1*dim2) - c1*dim2
    ! coordinates through a cartesian Virtual toplogy based on new_comm
    dims(0) = dim0; dims(1) = dim1;  periods(0)=.TRUE.; periods(1)=.TRUE.
    CALL MPI_CART_CREATE(new_comm, 2, dims, periods, .FALSE., comm_cart, ierror) ! No reorder !!!
    CALL MPI_CART_COORDS(comm_cart, my_new_rank, 2, coords, ierror)
    ! coparison of the results
    IF (c0 .NE. coords(0)) WRITE(*,'(''NEWrank='',I3,'', ERROR in coords(0): '',I2,'' != '',I2)') my_new_rank, c0, coords(0)
    IF (c1 .NE. coords(1)) WRITE(*,'(''NEWrank='',I3,'', ERROR in coords(0): '',I2,'' != '',I2)') my_new_rank, c1, coords(1)

! Ring in direction 0, i.e., with different c0 (and same other coord(s))
! ----------------------------------------------------------------------

    ! Source code without Cartesian Virtual Toplogy
    left__coord0 = mod(c0-1+dim0, dim0);   left__rank = left__coord0*dim1+c1
    right_coord0 = mod(c0+1,      dim0);   right_rank = right_coord0*dim1+c1
    ! coordinates through a cartesian Virtual toplogy based on new_comm
    ! right=..., left=... should be substituted by one call to MPI_Cart_shift():
    CALL MPI_CART_SHIFT(comm_cart, 0, 1, left, right, ierror)  ! dir=0 !!!
    ! coparison of the results
    IF (left__rank.NE.left ) WRITE(*,'(''DIR=0, NEWrank='',I3,'', ERROR in left:  '',I2,'' != '',I2)') my_new_rank,left__rank,left
    IF (right_rank.NE.right) WRITE(*,'(''DIR=0, NEWrank='',I3,'', ERROR in right: '',I2,'' != '',I2)') my_new_rank,right_rank,right

    sum = 0
    snd_buf = my_new_rank
  
    DO i=1, dim0
      ! without Cartesian Virtuzal Topology, comm_cart must be substituted by new-comm
      CALL MPI_ISSEND(snd_buf, 1, MPI_INT, right, to_right, comm_cart, request, ierror)
      CALL MPI_RECV(rcv_buf, 1, MPI_INT, left, to_right, comm_cart, status, ierror)
      CALL MPI_WAIT(request, status, ierror)
      snd_buf = rcv_buf
      sum = sum + rcv_buf
    END DO
    WRITE(*,'(''DIR=0, RANK world: '',I3,'' new: '',I3,'' -- coords(0): '',I2,'' (1): '',I2,'// &
     &                                            ''' -- left= '',I3,'' right= '',I3,'' -- sum= '',I4)')& 
     &               my_world_rank, my_new_rank,coords(0),coords(1),   left,right,        sum

! Ring in direction 1, i.e., with different c1 (and same other coord(s))
! ----------------------------------------------------------------------

    ! Source code without Cartesian Virtual Toplogy
    left__coord1 = mod(c1-1+dim1, dim1);   left__rank = c0*dim1 + left__coord1
    right_coord1 = mod(c1+1,      dim1);   right_rank = c0*dim1 + right_coord1
    ! coordinates through a cartesian Virtual toplogy based on new_comm
    ! right=..., left=... should be substituted by one call to MPI_Cart_shift():
    CALL MPI_CART_SHIFT(comm_cart, 1, 1, left, right, ierror)  ! dir=1 !!!
    ! coparison of the results
    IF (left__rank.NE.left ) WRITE(*,'(''DIR=1, NEWrank='',I3,'', ERROR in left:  '',I2,'' != '',I2)') my_new_rank,left__rank,left
    IF (right_rank.NE.right) WRITE(*,'(''DIR=1, NEWrank='',I3,'', ERROR in right: '',I2,'' != '',I2)') my_new_rank,right_rank,right

    sum = 0
    snd_buf = my_new_rank
  
    DO i=1, dim1
      ! without Cartesian Virtuzal Topology, comm_cart must be substituted by new-comm
      CALL MPI_ISSEND(snd_buf, 1, MPI_INT, right, to_right, comm_cart, request, ierror)
      CALL MPI_RECV(rcv_buf, 1, MPI_INT, left, to_right, comm_cart, status, ierror)
      CALL MPI_WAIT(request, status, ierror)
      snd_buf = rcv_buf
      sum = sum + rcv_buf
    END DO
    WRITE(*,'(''DIR=1, RANK world: '',I3,'' new: '',I3,'' -- coords(0): '',I2,'' (1): '',I2,'// &
     &                                            ''' -- left= '',I3,'' right= '',I3,'' -- sum= '',I4)')& 
     &               my_world_rank, my_new_rank,coords(0),coords(1),   left,right,        sum
  END IF

  DEALLOCATE(ranks) ! can be done already directly after the call to MPI_Group_incl()

  CALL MPI_FINALIZE(ierror)

END PROGRAM
