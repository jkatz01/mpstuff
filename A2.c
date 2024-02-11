#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define RAND_INC 6
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

void print_array(int* data, int size, const char* prefix) {
	printf("%s", prefix);
	for (int i = 0; i < size; i++) {
		printf("%3d", data[i]);
		if (i < size-1) {
			printf(", ");
		}
	}
	printf("\n");
}

int find_b_index(int max_a, int* data_b, int start, int b_size) {
	// uses binary search to find the greated j
	// that is smaller or equal to max_a
	printf("Max: %d Start: %d Size: %d", max_a, start, b_size);
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

	while (low < b_size) {
		if (data_b[low] <= max_a) {
			printf("\neq\n");
			low++;
		}
		else {
			break;
		}
	}

	// this returns -1 if there is no element
	// in B that is <= max_a
	printf("Low: %d\n", low-1);
	return low-1;
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
	int*	data_b_displs;

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


	data_a_sizes   = (int*)malloc(num_procs * sizeof(int));
	data_a_indices = (int*)malloc(num_procs * sizeof(int));
	data_b_sizes   = (int*)malloc(num_procs * sizeof(int));
	data_b_indices = (int*)malloc(num_procs * sizeof(int));
	data_b_displs  = (int*)malloc(num_procs * sizeof(int));
	
	if (my_rank == FIRST) {

		data_a = (int*)malloc(data_size * sizeof(int));
		data_b = (int*)malloc(data_size * sizeof(int));

		generate_arrays(data_a, data_b, data_size);

		partition_array(data_size, num_procs, data_a_sizes);

		int sum = 0;
		int i = 0;
		for (i = 0; i < num_procs; i++) {
			// calculate A indices
			data_a_indices[i] = sum;
			sum += data_a_sizes[i];
		}
		// Insert given message to appropriate location
		// Change to MPI Broadcast?
		//
		print_array(data_a, data_size, "A: ");
		print_array(data_b, data_size, "B: ");

		int previous_index = -1; //j()
		for (i = 0; i < num_procs; i++) {
			int cur_idx = data_a_indices[i];
			int cur_size = data_a_sizes[i];
			data_b_displs[i] = previous_index + 1; //displacements should start at 0 but indices are the max indices
			data_b_indices[i] = find_b_index(data_a[cur_idx + cur_size - 1], data_b, previous_index + 1, data_size);
			previous_index = data_b_indices[i];
		}

		int prev = -1; //first size includes 0
		for (i = 0; i < num_procs; i++) {
			data_b_sizes[i] = data_b_indices[i] - prev;
			prev = data_b_indices[i];
		}

		
		print_array(data_b_indices, num_procs, "j() indices: ");
		print_array(data_b_sizes, num_procs, "j() sizes: ");

	}
	
	MPI_Bcast(data_a_sizes, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_b_sizes, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	int my_a_size = data_a_sizes[my_rank];
	int my_b_size = data_b_sizes[my_rank];
	int* recv_a_buffer = (int*)malloc(my_a_size);
	int* recv_b_buffer = (int*)malloc(my_b_size);

	MPI_Scatterv(data_a, data_a_sizes, data_a_indices, MPI_INT, recv_a_buffer, my_a_size, MPI_INT, FIRST, MPI_COMM_WORLD);
	MPI_Scatterv(data_b, data_b_sizes, data_b_displs, MPI_INT, recv_b_buffer, my_b_size, MPI_INT, FIRST, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	printf("[%d]: ", my_rank);
	print_array(recv_a_buffer, my_a_size, "My A: ");
	printf("[%d]: ", my_rank);
	print_array(recv_b_buffer, my_b_size, "My B: ");

	
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
