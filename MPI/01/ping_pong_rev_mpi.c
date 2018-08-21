#include<stdio.h>
#include<mpi.h>

int main(){
        int n;
		int i;
        int my_rank, num_procs;

        MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		
		double rank0 = 123.0, rank1;
		MPI_Status stat;
		
		if (my_rank == 0){
			MPI_Send(&rank0, 1, MPI_INT, 1, 17, MPI_COMM_WORLD);
			MPI_Recv(&rank1, 1, MPI_INT, 1, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		}
		else{
			MPI_Recv(&rank1, 1, MPI_INT, 0, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&rank1, 1, MPI_INT, 0, 23, MPI_COMM_WORLD);
		}
			  
		double start = MPI_Wtime();
		
		for(i = 0; i < 50; ++i){
				if (my_rank == 0){
					MPI_Ssend(&rank0, 1, MPI_INT, 1, 17, MPI_COMM_WORLD);
					MPI_Recv(&rank1, 1, MPI_INT, 1, 23, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				}
				else{
					MPI_Recv(&rank1, 1, MPI_INT, 0, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Ssend(&rank1, 1, MPI_INT, 0, 23, MPI_COMM_WORLD);
				}
			  
		}
		if (my_rank == 0) {
			double end = MPI_Wtime();
			double time = end - start;
			double latency = time/(2.0*50)*1e6;
			printf("Calculate time %lf us.\n", latency);
			printf("The bandwidth is %f Bytes/us.\n", 8.0/latency);
		}
		
		MPI_Finalize();
        return 0;
}
