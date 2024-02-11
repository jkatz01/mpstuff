#include <stdio.h>
#include "mpi.h"

#define SIZE_D 6
#define SIZE_C (SIZE_A + SIZE_B)

int main (int argc, char *argv[]) {
	int	my_rank;
	int	p;
	int	source;
	int	dest;
	int	tag = 0;
	
	// Randomly generate these arrays
	// Can we do this with shared memory?
	int	data_a[SIZE_D] = {12, 13, 26, 40, 55, 71};
	int	data_b[SIZE_D] = { 5, 19, 21, 50, 51, 90};
	int	sum = 0;

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

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


	if (my_rank == 0) {
		// Insert given message to appropriate location
		for (source = 1; source < p; source++) {
			MPI_Recv(&sum, SIZE_D, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			printf("Sum: %d\n", sum);
		}
	}
	else {
		for (int i = 0; i < SIZE_D; i++) {
			sum += data_a[i];
		}
		dest = 0;
		MPI_Send(&sum, SIZE_D, MPI_INT, dest, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;

}

int binary_search(const int *data, int start, int end) {
	int index = 0;
	return index;
}
