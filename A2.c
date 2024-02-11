#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define RAND_INC 7
#define FIRST    0

enum Tags {
	TAG_CHUNK_START,
	TAG_CHUNK_END
};

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

int find_b_index(int max_a, int* data_b, int start, int b_size) {
	// uses binary search to find the greated j
	// that is smaller or equal to max_a
	int middle, low, high;
	low = start; // start is the previous j()
	high = b_size - 1;
	while (low < high) {
		middle = low + (high - low) / 2;

		if (data_b[middle] == max_a) {
			return middle;
		}
		if (data_b[middle] < max_a + 1) {
			low = middle + 1;
		}
		else {
			high = middle;
		}
	}
	// this returns -1 if there is no element
	// in B that is <= max_a
	return low - 1;
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
	int*	data_a_indices;
	int*	data_b_indices; //equivalent to j() in algorithm

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
	data_a_indices = (int*)malloc(num_procs * sizeof(int));
	data_b_indices = (int*)malloc(num_procs * sizeof(int));
	
	if (my_rank == FIRST) {

		data_a = (int*)malloc(data_size * sizeof(int));
		data_b = (int*)malloc(data_size * sizeof(int));

		generate_arrays(data_a, data_b, data_size);

		partition_array(data_size, num_procs, data_a_sizes);

		int sum = 0;
		for (int i = 0; i < num_procs; i++) {
			// calculate A indices
			data_a_indices[i] = sum;
			sum += data_a_sizes[i];
		}
		// Insert given message to appropriate location
		// Change to MPI Broadcast?
		//
		print_array(data_a, data_size);
		print_array(data_b, data_size);

		for (int j = 0; j < num_procs; j++) {
			int cur_idx = data_a_indices[j];
			int cur_size = data_a_sizes[j];
			data_b_indices[j] = find_b_index(data_a[cur_idx + cur_size], data_b, cur_idx, data_size);
		}

		print_array(data_b_indices, num_procs);

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(data_a_sizes, data_size, MPI_INT, 0, MPI_COMM_WORLD);
	
	int my_a_size = data_a_sizes[my_rank];
	int* recv_a_buffer = (int*)malloc(my_a_size);

	MPI_Scatterv(data_a, data_a_sizes, data_a_indices, MPI_INT, recv_a_buffer, my_a_size, MPI_INT, FIRST, MPI_COMM_WORLD);

	printf("[%d]: ", my_rank);
	print_array(recv_a_buffer, my_a_size);

	
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
