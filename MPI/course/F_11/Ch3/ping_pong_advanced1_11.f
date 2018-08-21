      PROGRAM pingpong

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                              C
C This file has been written as a sample solution to an        C 
C exercise in a course given at the High Performance           C
C Computing Centre Stuttgart (HLRS).                           C
C The examples are based on the examples in the MPI course of  C
C the Edinburgh Parallel Computing Centre (EPCC).              C
C It is made freely available with the understanding that      C 
C every copy of this file must include this header and that    C
C HLRS and EPCC take no responsibility for the use of the      C
C enclosed teaching material.                                  C
C                                                              C
C Authors: Joel Malard, Alan Simpson,            (EPCC)        C
C          Rolf Rabenseifner, Traugott Streicher (HLRS)        C
C                                                              C
C Contact: rabenseifner@hlrs.de                                C 
C                                                              C  
C Purpose: A program to try MPI_Ssend and MPI_Recv.            C
C                                                              C
C Contents: F-Source                                           C
C                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT NONE

      INCLUDE "mpif.h"

      INTEGER proc_a
      PARAMETER(proc_a=0)
                
      INTEGER proc_b
      PARAMETER(proc_b=1)                

      INTEGER ping
      PARAMETER(ping=17)
        
      INTEGER pong
      PARAMETER(pong=23)        

      INTEGER number_of_messages 
      PARAMETER (number_of_messages=50)

      INTEGER length
      PARAMETER (length=1)
 
      DOUBLE PRECISION start, finish, time
      INTEGER status(MPI_STATUS_SIZE)
   
      REAL buffer(length)

      INTEGER i

      INTEGER ierror, my_rank, size


      CALL MPI_INIT(ierror)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)

      IF (my_rank .EQ. proc_a) THEN
            CALL MPI_SEND(buffer, length, MPI_REAL,
     &                     proc_b, ping, MPI_COMM_WORLD, ierror)
            CALL MPI_RECV(buffer, length, MPI_REAL,
     &                    proc_b, pong, MPI_COMM_WORLD,
     &                    status, ierror)
      ELSE IF (my_rank .EQ. proc_b) THEN
            CALL MPI_RECV(buffer, length, MPI_REAL,
     &                    proc_a, ping, MPI_COMM_WORLD,
     &                    status, ierror)
            CALL MPI_SEND(buffer, length, MPI_REAL,
     &                     proc_a, pong, MPI_COMM_WORLD, ierror)
      END IF

      start = MPI_WTIME()

      DO i = 1, number_of_messages

         IF (my_rank .EQ. proc_a) THEN
            CALL MPI_SEND(buffer, length, MPI_REAL,
     &                     proc_b, ping, MPI_COMM_WORLD, ierror)
            CALL MPI_RECV(buffer, length, MPI_REAL,
     &                    proc_b, pong, MPI_COMM_WORLD,
     &                    status, ierror)
         ELSE IF (my_rank .EQ. proc_b) THEN
            CALL MPI_RECV(buffer, length, MPI_REAL,
     &                    proc_a, ping, MPI_COMM_WORLD,
     &                    status, ierror)
            CALL MPI_SEND(buffer, length, MPI_REAL,
     &                     proc_a, pong, MPI_COMM_WORLD, ierror)
         END IF

      END DO

      finish = MPI_WTIME()

      IF (my_rank .EQ. proc_a) THEN

         time = finish - start

         WRITE(*,*) 'Time for one message:', 
     &           time/(2*number_of_messages)*1e6, ' micro seconds'

      END IF

      CALL MPI_FINALIZE(ierror)

      END
