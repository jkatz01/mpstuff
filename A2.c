#include <stdio.h>
#include "mpi.h"

#define SIZE_A 6
#define SIZE_B 6
#define SIZE_C (SIZE_A + SIZE_B)

int main (int argc, char *argv[]) {
	int	my_rank;
	int	p;
	int	source;
	int	dest;
	int	tag = 0;

	int	data_a[SIZE_A] = {12, 56, 72, 13, 26, 40};
	int	data_b[SIZE_B] = {87, 21, 19, 50, 51, 90};
	int	sum = 0;

	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (my_rank == 0) {
		for (source = 1; source < p; source++) {
			MPI_Recv(&sum, SIZE_A, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			printf("Sum: %d\n", sum);
		}
	}
	else {
		for (int i = 0; i < SIZE_A; i++) {
			sum += data_a[i];
		}
		dest = 0;
		MPI_Send(&sum, SIZE_A, MPI_INT, dest, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;

}
