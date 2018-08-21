#include<stdio.h>
#include<mpi.h>

int main(){
        int n;
		int i;
        int my_rank, num_procs;

        MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		
		int rank0=1, rank1=0;
		
		
		MPI_Status stat;
		MPI_Request request;
		int tag = 17;
		/*int sum[num_procs];*/
		int sum = 0;
		
		for (i = 0; i< num_procs; ++i){
			int destination = (my_rank + 1) % num_procs;
			int source = (my_rank - 1 + num_procs) % num_procs;
			MPI_Issend(&rank0, 1, MPI_INT, destination, tag, MPI_COMM_WORLD, &request);
			MPI_Recv(&rank1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &stat);
			
			/*MPI_Ssend(&rank0, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
			MPI_Irecv(&rank1, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &request);*/
			
			MPI_Wait(&request, &stat);
			rank0 = rank1;
			sum += rank1;
			
			/*printf("The sum is %i.\n", sum);*/
		}
		printf("The sum is %i.\n", sum);
		/*for (i = 0; i< num_procs; ++i){
			printf("The sum is %i.\n", sum[i]);
		}*/
		MPI_Finalize();
        return 0;
}