#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define RAND_INC 7

enum Tags {
	TAG_CHUNK_START,
	TAG_CHUNK_END
};

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

/* chunk_sizes must have total_procs - 1 elements because 
 * process 0 is not being used				*/
void partition_array(int data_size, int num_procs, int* chunks) {
	/* chunk_sizes must have total num of procs elements
	
	*  determine the size of each 
	*  chunk, and then incrementally
	*  add to them until the remainder
	*  is 0					*/
	
	int remainder = data_size % num_procs;
	int initial = data_size / num_procs;

	for (int i = 0; i < num_procs; i++) {
		chunks[i] = initial;
		if (remainder > 0) {
			chunks[i] += 1;
			remainder--;
		}
	}

}

void print_array(int* data, int size) {
	for (int i = 0; i < size; i++) {
		printf("%d", data[i]);
		if (i < size-1) {
			printf(", ");
		}
	}
	printf("\n");
}

int main (int argc, char *argv[]) {
	int	my_rank;
	int	num_procs;
	int	source;
	int	dest;
	int	data_size;
	int	k;
	int*	data_a;
	int*	data_b;
	int*	data_a_sizes;
	int*	data_b_sizes;

	if (argc != 2) {
		data_size = 10;
	}
	else {
		data_size = atoi(argv[1]);
	}

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


	data_a_sizes = (int*)malloc(num_procs * sizeof(int));
	
	if (my_rank == 0) {

		data_a = (int*)malloc(data_size * sizeof(int));
		data_b = (int*)malloc(data_size * sizeof(int));

		generate_arrays(data_a, data_b, data_size);

		partition_array(data_size, num_procs, data_a_sizes);

		// Insert given message to appropriate location
		// Change to MPI Broadcast?

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(data_a_sizes, data_size, MPI_INT, 0, MPI_COMM_WORLD);
	print_array(data_a_sizes, num_procs);

	MPI_Finalize();

	return 0;

	// TODO: use proper logging instead of printing to console

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
}
