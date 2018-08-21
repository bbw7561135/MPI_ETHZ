PROGRAM ibarrier

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
! Authors: Rolf Rabenseifner (HLRS)                            !
!                                                              !
! Contact: rabenseifner@hlrs.de                                !
!                                                              !
! Purpose: A program to try MPI_Ibarrier.                      !
!                                                              !
! Contents: C-Source                                           !
!                                                              !
!**************************************************************!


  USE mpi_f08

  IMPLICIT NONE

  INTEGER :: my_rank, size
! ...in the role as sending process
  INTEGER :: snd_buf_A, snd_buf_B, snd_buf_C, snd_buf_D
  INTEGER :: dest, number_of_dests=0
  LOGICAL :: snd_finished=.FALSE.
  TYPE(MPI_Request) :: snd_rq(0:3)
  INTEGER :: total_number_of_dests ! only for verification, should be removed in real applications 
                                   ! Caution: total_number_of_dests may be less than 4, see IF-statements below
! ...in the role as receiving process
  INTEGER :: rcv_buf
  TYPE(MPI_Request) :: ib_rq
  LOGICAL :: ib_finished=.FALSE., rcv_flag
  TYPE(MPI_Status) :: rcv_sts
  INTEGER :: number_of_recvs=0, total_number_of_recvs ! only for verification, should be removed in real applications

  INTEGER :: round=0 ! only for verification, should be removed in real applications

  CALL MPI_Init()
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size)

! ...in the role as sending process
  dest = my_rank+1
  IF ((dest>=0) .AND. (dest<size)) THEN
    snd_buf_A = 1000*my_rank + dest   ! must not be modified until send-completion with TEST or WAIT 
    CALL MPI_Issend(snd_buf_A,1,MPI_INTEGER, dest,222,MPI_COMM_WORLD, snd_rq(number_of_dests))
    write(*,'(''A rank: '',I3,'' - sending_: message '',I6.6,'' from '',I3,'' to '',I3)') my_rank, snd_buf_A, my_rank, dest
    number_of_dests = number_of_dests + 1
  END IF
  dest = my_rank-2
  IF ((dest>=0) .AND. (dest<size)) THEN
    snd_buf_B = 1000*my_rank + dest   ! must not be modified until send-completion with TEST or WAIT 
    CALL MPI_Issend(snd_buf_B,1,MPI_INTEGER, dest,222,MPI_COMM_WORLD, snd_rq(number_of_dests))
    write(*,'(''A rank: '',I3,'' - sending_: message '',I6.6,'' from '',I3,'' to '',I3)') my_rank, snd_buf_B, my_rank, dest
    number_of_dests = number_of_dests + 1
  END IF
  dest = my_rank+4
  IF ((dest>=0) .AND. (dest<size)) THEN
    snd_buf_C = 1000*my_rank + dest   ! must not be modified until send-completion with TEST or WAIT 
    CALL MPI_Issend(snd_buf_C,1,MPI_INTEGER, dest,222,MPI_COMM_WORLD, snd_rq(number_of_dests))
    write(*,'(''A rank: '',I3,'' - sending_: message '',I6.6,'' from '',I3,'' to '',I3)') my_rank, snd_buf_C, my_rank, dest
    number_of_dests = number_of_dests + 1
  END IF
  dest = my_rank-7
  IF ((dest>=0) .AND. (dest<size)) THEN
    snd_buf_D = 1000*my_rank + dest   ! must not be modified until send-completion with TEST or WAIT 
    CALL MPI_Issend(snd_buf_D,1,MPI_INTEGER, dest,222,MPI_COMM_WORLD, snd_rq(number_of_dests))
    write(*,'(''A rank: '',I3,'' - sending_: message '',I6.6,'' from '',I3,'' to '',I3)') my_rank, snd_buf_D, my_rank, dest
    number_of_dests = number_of_dests + 1
  END IF
  DO WHILE(.NOT. ib_finished)
!   ...in the role as receiving process
!       MPI_IPROBE(MPI_ANY_SOURCE); If there is a message then MPI_RECV for this one message:
        CALL MPI_Iprobe(______________,___,______________, ________, _________________) 
        IF(rcv_flag) THEN
          CALL MPI_Recv(rcv_buf,1,MPI_INTEGER, MPI_ANY_SOURCE,222,MPI_COMM_WORLD, rcv_sts)
          write(*,'(''A rank: '',I3,'' - received: message '',I6.6,'' from '',I3,'' to '',I3)') &
&                            my_rank,                     rcv_buf, rcv_sts%MPI_SOURCE, my_rank
          number_of_recvs = number_of_recvs + 1  ! only for verification, should be removed in real applications
        END IF
!   ...in the role as sending process:
!      The following lines make only sense as long as not all MPI_ISSENDs are finished.
    IF(.NOT. snd_finished) THEN
      ! Check whether all MPI_ISSENDs are finished
      CALL MPI_Te_____(_______________, ______, ____________, ___________________)
      ! if all MPI_ISSENDs are finished then call MPI_IBARRIER
      IF(____________) THEN   ! ...i.e., the first time, i.e., only once
        CALL MPI_________(______________, _____)
      END IF
    END IF
!   ...loop until MPI_IBARRIER finished (i.e. all processes signaled that all receives are called)
    IF(____________) THEN ! the test whether the MPI_IBARRIER is finished
                          ! can be done only if MPI_IBARRIER is already started.
                          ! This ist true as soon snd_finished is true
      CALL MPI_____(_____, ___________, MPI_STATUS_IGNORE)
    END IF
  END DO

  ! only for verification, should be removed in real applications:
  CALL MPI_Reduce(number_of_dests, total_number_of_dests, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD)
  CALL MPI_Reduce(number_of_recvs, total_number_of_recvs, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD)
  IF (my_rank .EQ. 0) THEN
    write(*,'(''B #sends= '',I5,''  /  #receives= '',I5)') total_number_of_dests, total_number_of_recvs
    IF (total_number_of_dests .NE. total_number_of_recvs) THEN
      write(*,'(''C ERROR !!!! Wrong number of receives'')')
    END IF
  END IF

  CALL MPI_Finalize()

END PROGRAM
