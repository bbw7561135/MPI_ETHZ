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

      INTEGER start_length 
      PARAMETER (start_length=8)

      INTEGER length_factor 
      PARAMETER (length_factor=64)

      INTEGER max_length                ! 2 Mega 
      PARAMETER (max_length=2097152)

      INTEGER number_package_sizes 
      PARAMETER (number_package_sizes=4)

      INTEGER i, j, size_of_real

      INTEGER length
 
      DOUBLE PRECISION start, finish, time, transfer_time
      INTEGER status(MPI_STATUS_SIZE)
   
      REAL buffer(max_length)

      INTEGER ierror, my_rank, size


      CALL MPI_INIT(ierror)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
      CALL MPI_TYPE_EXTENT(MPI_REAL, size_of_real, ierror) 

      IF (my_rank .EQ. proc_a) THEN
         WRITE (*,*) "message size   transfertime    bandwidth"
      END IF

      length = start_length

      DO j = 1, number_package_sizes

         IF (my_rank .EQ. proc_a) THEN
               CALL MPI_SSEND(buffer, length, MPI_REAL,
     &                        proc_b, ping, MPI_COMM_WORLD, ierror)
               CALL MPI_RECV(buffer, length, MPI_REAL,
     &                       proc_b, pong, MPI_COMM_WORLD,
     &                       status, ierror)
         ELSE IF (my_rank .EQ. proc_b) THEN
               CALL MPI_RECV(buffer, length, MPI_REAL,
     &                       proc_a, ping, MPI_COMM_WORLD,
     &                       status, ierror)
               CALL MPI_SSEND(buffer, length, MPI_REAL,
     &                        proc_a, pong, MPI_COMM_WORLD, ierror)
         END IF
     
         start = MPI_WTIME()
     
         DO i = 1, number_of_messages
     
            IF (my_rank .EQ. proc_a) THEN
               CALL MPI_SSEND(buffer, length, MPI_REAL,
     &                        proc_b, ping, MPI_COMM_WORLD, ierror)
               CALL MPI_RECV(buffer, length, MPI_REAL,
     &                       proc_b, pong, MPI_COMM_WORLD,
     &                       status, ierror)
            ELSE IF (my_rank .EQ. proc_b) THEN
               CALL MPI_RECV(buffer, length, MPI_REAL,
     &                       proc_a, ping, MPI_COMM_WORLD,
     &                       status, ierror)
               CALL MPI_SSEND(buffer, length, MPI_REAL,
     &                        proc_a, pong, MPI_COMM_WORLD, ierror)
            END IF
     
         END DO
     
         finish = MPI_WTIME()
     
         IF (my_rank .EQ. proc_a) THEN
     
            time = finish - start
            transfer_time = time / (2 * number_of_messages)
     
            WRITE(*,*) length*size_of_real, 'bytes  ',
     &                    transfer_time*1e6, 'usec  ',
     &                    1e-6*length*size_of_real/transfer_time, 'MB/s'
     
         END IF

         length = length * length_factor

      END DO

      CALL MPI_FINALIZE(ierror)

      END
