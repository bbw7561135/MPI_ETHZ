#include<stdio.h>
#include<mpi.h>

int main(){
        int n;
	int i;
        int my_rank, num_procs;

        MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		int rank0 = 123, rank1;
		MPI_Status stat;
		
		double start = MPI_Wtime();
		for(i = 0; i < 50; ++i){
			if (i == 49 && my_rank == 0){
					double end = MPI_Wtime();
					double time = end - start;
					printf("Calculate time %f us\n", time/(2.0*50)*1e6);	
			}
			else{
				if (my_rank == 0){
					MPI_Send(&rank0, 1, MPI_INT, 1, 17, MPI_COMM_WORLD);
					MPI_Recv(&rank1, 1, MPI_INT, 1, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				}
				else{
					MPI_Recv(&rank1, 1, MPI_INT, 0, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send(&rank1, 1, MPI_INT, 0, 23, MPI_COMM_WORLD);
				}
			}
			  
		}
		
		
		MPI_Finalize();
        return 0;
}
