#include<stdio.h>
#include<mpi.h>

int main(){
        int n;
        int my_rank, num_procs;

        MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        if (my_rank == 0){
			printf("Hello, World!\n");  
			printf("I am %i of %i\n", my_rank, num_procs);
		}
		else{
			printf("I am %i of %i\n", my_rank, num_procs);
		}
	MPI_Finalize();  
        return 0;
}
