        PROGRAM ring

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C This file has been written as a sample solution to an exercise in a 
C course given at the Edinburgh Parallel Computing Centre. It is made
C freely available with the understanding that every copy of this file
C must include this header and that EPCC takes no responsibility for
C the use of the enclosed teaching material.
C
C Authors:    Alan Simpson, Joel Malard, Rolf Rabenseifner
C
C Contact:    epcc-tec@epcc.ed.ac.uk, rabenseifner@hlrs.de 
C
C Purpose:    A program to test non-blocking point-to-point
C             communications.
C
C Contents:   F source code.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        IMPLICIT NONE

	INCLUDE "mpif.h"

	INTEGER tag_to_right
	PARAMETER(tag_to_right=201)

        INTEGER ierror, snd_buf, my_rank, size
	INTEGER right, left
	INTEGER rcv_buf, sum, i

	INTEGER MAXDIMS
	PARAMETER (MAXDIMS=2) 

	INTEGER new_comm
	INTEGER dims(MAXDIMS)
	LOGICAL periods(MAXDIMS), reorder
	INTEGER coords(MAXDIMS)

	INTEGER recv_status(MPI_STATUS_SIZE)

        CALL MPI_INIT(ierror)

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

C SET CARTESIAN TOPOLOGY.
 
	dims(1) = 0
	dims(2) = 0
	CALL MPI_DIMS_CREATE(size, MAXDIMS, dims, ierror) 
	periods(1) = .TRUE.
	periods(2) = .FALSE.
	reorder = .TRUE.

	CALL MPI_CART_CREATE(MPI_COMM_WORLD, MAXDIMS, dims, periods,
     &                       reorder, new_comm, ierror)
        CALL MPI_COMM_RANK(new_comm, my_rank, ierror)
	CALL MPI_CART_COORDS(new_comm,my_rank, MAXDIMS,coords, ierror) 

C GET NEAREST NEIGHBOUR RANKS.

	CALL MPI_CART_SHIFT(new_comm, 0, 1, left, right, ierror)

C COMPUTE SUM 

	sum = 0
	snd_buf = my_rank

	DO I=1, DIMS(1) 

          CALL MPI_SENDRECV(snd_buf, 1,MPI_INTEGER, right, tag_to_right,
     &                      rcv_buf, 1,MPI_INTEGER, left,  tag_to_right,
     &                      new_comm, recv_status, ierror)

	  sum = sum + rcv_buf
	  snd_buf = rcv_buf

	END DO 

	WRITE(*,*) 'PE=',my_rank,' x=',coords(1),' y=',coords(2),
     &             '  sum=',sum

        CALL MPI_FINALIZE(ierror)

        STOP
        END
