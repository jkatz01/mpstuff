#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define RAND_INC 6
#define FIRST    0

/* Generates two arrays of semi-random numbers */
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

/* Partitions an array into as balanced as possible chunks */
void partition_array(int data_size, int num_procs, int* chunks) {
	int remainder	= data_size % num_procs;
	int initial	= data_size / num_procs;

	for (int i = 0; i < num_procs; i++) {
		chunks[i] = initial;
		if (remainder > 0) {
			chunks[i] += 1;
			remainder--;
		}
	}

}

/* Prints an array with a prefix */
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

/* Finds an index j that is the greatest value in data_b
 * that is >= max_a					*/
int find_b_index(int max_a, int* data_b, int start, int b_size) {
	// uses binary search to find the greated j
	// that is smaller or equal to max_a
	int middle, low, high;
	low = start; // start is the previous j()
	high = b_size - 1;
	while (low < high) {
		middle = low + (high - low) / 2;

		if (data_b[middle] == max_a) {
			break;
		}
		if (data_b[middle] < max_a + 1) {
			low = middle + 1;
		}
		else {
			high = middle;
		}
	}

	// Cover duplicates
	while (low < b_size) {
		if (data_b[low] <= max_a) {
			low++;
		}
		else {
			break;
		}
	}

	// this returns -1 if there is no element
	// in B that is <= max_a
	return low-1;
}

/* Merges two sorted arrays into one */
void merge_arrays(const int* data_a, const int* data_b, int* data_c, int a_size, int b_size) {
	int k = 0;
	int i = 0;
	int j = 0;

	while (i < a_size && j < b_size) {
		if (data_a[i] < data_b[j]) {
			data_c[k] = data_a[i];
			k++;
			i++;
		}
		else {
			data_c[k] = data_b[j];
			k++;
			j++;
		}
	}
	if (i == a_size) {
		while (j < b_size) {
			data_c[k] = data_b[j];
			k++;
			j++;
		}
	}
	if (j == b_size) {
		while (i < a_size) {
			data_c[k] = data_a[i];
			k++;
			i++;
		}
	}
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
	int*	data_c;
	int*	data_a_sizes;
	int*	data_b_sizes;
	int*	data_c_sizes;
	int*	data_a_indices;
	int*	data_b_indices;
	int*	data_b_displs;
	int*	data_c_displs;

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

	// This really should just be a struct lmao
	data_a_sizes   = (int*)malloc(num_procs * sizeof(int));
	data_a_indices = (int*)malloc(num_procs * sizeof(int));
	data_b_sizes   = (int*)malloc(num_procs * sizeof(int));
	data_b_indices = (int*)malloc(num_procs * sizeof(int));
	data_b_displs  = (int*)malloc(num_procs * sizeof(int));
	data_c_sizes   = (int*)malloc(num_procs * sizeof(int));
	data_c_displs  = (int*)malloc(num_procs * sizeof(int));
	
	if (my_rank == FIRST) {

		data_a = (int*)malloc(data_size * sizeof(int));
		data_b = (int*)malloc(data_size * sizeof(int));

		generate_arrays(data_a, data_b, data_size);

		partition_array(data_size, num_procs, data_a_sizes);

		int sum = 0;
		int i = 0;
		for (i = 0; i < num_procs; i++) {
			data_a_indices[i] = sum;
			sum		  += data_a_sizes[i];
		}

		int previous_index = -1; //j()
		for (i = 0; i < num_procs; i++) {
			int cur_idx		= data_a_indices[i];
			int cur_size		= data_a_sizes[i];
			data_b_displs[i]	= previous_index + 1; //displacements should start at 0 but indices are the max indices
			data_b_indices[i]	= find_b_index(data_a[cur_idx + cur_size - 1], data_b, previous_index + 1, data_size);
			previous_index		= data_b_indices[i];
		}

		int prev = -1; //first size includes 0
		for (i = 0; i < num_procs; i++) {
			data_b_sizes[i] = data_b_indices[i] - prev;
			prev		= data_b_indices[i];
		}
		
		// data c sizes
		int previous_c_displs = 0;
		for (i = 0; i < num_procs; i++) {
			int cur_size		= data_a_sizes[i] + data_b_sizes[i];
			data_c_sizes[i]		= cur_size;
			data_c_displs[i]	= previous_c_displs;
			previous_c_displs	+= cur_size;
		}

	}

	MPI_Bcast(data_a_sizes,   num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_b_sizes,   num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_c_sizes,   num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_a_indices, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_b_indices, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_b_displs,  num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(data_c_displs,  num_procs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);	

	
	int my_a_size		= data_a_sizes[my_rank];
	int my_b_size		= data_b_sizes[my_rank];
	int my_c_size		= my_a_size + my_b_size;
	int* recv_a_buffer	= (int*)malloc(my_a_size * sizeof(int));
	int* recv_b_buffer	= (int*)malloc(my_b_size * sizeof(int));

	MPI_Scatterv(data_a, data_a_sizes, data_a_indices, MPI_INT, recv_a_buffer, my_a_size, MPI_INT, FIRST, MPI_COMM_WORLD);
	MPI_Scatterv(data_b, data_b_sizes, data_b_displs,  MPI_INT, recv_b_buffer, my_b_size, MPI_INT, FIRST, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	int* my_data_c = (int*)malloc(my_c_size * sizeof(int));

	// Here's the actual parallel part of this problem
	// I don't know how to calculate j() in parallel without
	// copying the entire B array for every process
	merge_arrays(recv_a_buffer, recv_b_buffer, my_data_c, my_a_size, my_b_size);

	if (my_rank == FIRST) {
		// only the first element needs a valid receive buffer
		data_c = (int*)malloc((data_size * 2) * sizeof(int)); // A and B are same length
	}
	else {
		data_c = NULL;
	}
	MPI_Gatherv(my_data_c, my_c_size, MPI_INT, data_c, data_c_sizes, data_c_displs, MPI_INT, FIRST, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank == FIRST) {
		//print_array(data_c, data_size*2, "C: ");
	}
	MPI_Finalize();

	return 0;

}
