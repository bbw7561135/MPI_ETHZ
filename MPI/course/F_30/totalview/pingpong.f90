!   Program for totalview exercise
!   by Thomas Boenisch, HLRS
!

program pingpong

   include 'mpif.h'


   integer array(1000,1000)
   
   character*55 message

   integer rank, ierror

   integer status(MPI_STATUS_SIZE)
   
   integer i_buf(100000)

   
   call mpi_init(ierror)

   call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)

   
   if (rank == 0) then

      message = 'hello, this is the MPI-Course at the HLRS, April 2000 '
         
      do i = 1, 53 

		 call mpi_send(message,i,MPI_CHAR,1,4,MPI_COMM_WORLD,ierror)

		 call mpi_recv(ibuf,0,MPI_INTEGER,1,5,MPI_COMM_WORLD,status,ierror)

      end do

      call mpi_send(array,1000000,MPI_INTEGER,1,6,MPI_COMM_WORLD,ierror)

      call mpi_recv(array,1000000,MPI_INTEGER,1,7,MPI_COMM_WORLD,status,ierror)


   else

	  do i = 1, 53 

         message = '     '
		 
		 call mpi_recv(message,i,MPI_CHAR,0,5,MPI_COMM_WORLD,status,ierror)

		 call mpi_send(ibuf,0,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierror)

		 write(6,*) message

      end do

      call mpi_recv(ibuf,1000000,MPI_INTEGER,0,6,MPI_COMM_WORLD,status,ierror)

      call mpi_send(ibuf,1000000,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierror)

   endif
   
   call mpi_finalize(ierror)

end
