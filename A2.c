#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define RAND_INC 7

int binary_search(const int *data, int start, int end);
void generate_arrays(int* array_a, int* array_b, int size);

int main (int argc, char *argv[]) {
	int	my_rank;
	int	procs;
	int	source;
	int	dest;
	int	tag = 0;
	int	data_size;
	
	// Randomly generate these arrays
	// Can we do this with shared memory?

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// Partition A into r groups
	// each with k = logn elements
	// Group 1: A[1]	... A[k]
	// Group 2: A[k+1]	... A[2k]
	// Group i: A[(i-1)k+1]	... A[ik]
	// Group r: A[(r-1)k+1]	... A[rk]
	//
	// Find r integers j(1) ... j(r)
	// j(1) is greatest index so A[k]  >= B[j(1)]
	// j(2) is greatest index so A[2k] >= B[j(2)]
	// j(i) is greatest index so A[ik] >= B[j(i)]
	// j(r) is greatest index so A[rk] >= B[j(r)]
	//
	// This partitions B into r groups
	// b[1]	      ...b[j(1)], b[j(1)+1]  ...b[j(2)]
	// b[j(i-1)+1]...b[j(i)], b[j(r-1)+1]...b[j(r)]
	//
	// Assign processor i to merge group i of A 
	//			     & group i of B
	//
	// This guarantees that the elments of B have
	// reached their final position in C(1:2n)

	int* data_a;
	int* data_b;
	data_size = 10;

	if (my_rank == 0) {

		data_a = (int*)malloc(data_size * sizeof(int));
		data_b = (int*)malloc(data_size * sizeof(int));

		generate_arrays(data_a, data_b, data_size);
		
		// Insert given message to appropriate location
		for (source = 1; source < procs; source++) {
			MPI_Send(&data_a, data_size, MPI_INT, source, tag, MPI_COMM_WORLD);
			//MPI_Recv(&sum, SIZE_D, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}
	else {
		//MPI_Send(&sum, SIZE_D, MPI_INT, dest, tag, MPI_COMM_WORLD);
		MPI_Recv(&data_a, data_size, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		for (int i = 0; i < data_size; i++) {
			printf("%d, ", data_a[i]);
		}
		printf("\n");
	}

	MPI_Finalize();

	return 0;

}

int binary_search(const int *data, int start, int end) {
	int index = 0;
	return index;
}

void generate_arrays(int* data_a, int* data_b, int size) {
	int cur_a = 0;
	int cur_b = 0;
	for (int i = 0; i < size; i++) {
		cur_a += rand() % RAND_INC;
		cur_b += rand() % RAND_INC;
		data_a[i] = cur_a;
		data_b[i] = cur_b;
	}
}
