/*   Program for totalview exercise
     by Thomas Boenisch, HLRS
*/


#include <mpi.h>

main(int argc,char *argv[])
{

   int array[1000][1000];

   char message[55] = "hello, this is the MPI-Course at the HLRS, April 2000";
         
   char message_buf[60];

   int rank, i;

   MPI_Status status;
   
   int ibuf;

   
   MPI_Init(&argc,&argv);

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   
   if (rank == 0){
     
     for(i=0;i<54;i++){
       
       MPI_Send(message,i,MPI_CHAR,1,4,MPI_COMM_WORLD);
	 
       MPI_Recv(&ibuf,0,MPI_INT,1,5,MPI_COMM_WORLD,&status);

     }

     MPI_Send(array,1000000,MPI_INT,1,6,MPI_COMM_WORLD);

     MPI_Recv(array,1000000,MPI_INT,1,7,MPI_COMM_WORLD,&status);

   }
   else{

     for(i=0;i<60;i++) message_buf[i]='\0';

     for(i=0;i<54;i++){
	
       MPI_Recv(message_buf,i,MPI_CHAR,0,5,MPI_COMM_WORLD,&status);

       MPI_Send(&ibuf,0,MPI_INT,0,4,MPI_COMM_WORLD);

       printf("%s\n", message_buf);

     }

     MPI_Recv(&ibuf,1000000,MPI_INT,0,6,MPI_COMM_WORLD,&status);
       
     MPI_Send(&ibuf,1000000,MPI_INT,0,7,MPI_COMM_WORLD);

   }
   
   MPI_Finalize();

 }


