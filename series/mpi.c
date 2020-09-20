#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#include "../funcs.h"

#define ROW_SIZE	480
#define	COLUMN_SIZE	480

#define REPEAT_TIMES 20
#define MAX_TIMES 500
#define TERMCHECK_TIMES 10

int main() {
	
	MPI_Comm	comm;
	int ndims, reorder, periods[2], dim_size[2];
	int rank,size;
	int * cells = NULL; 
	int * np_cells = NULL;
	int * temp = NULL;
	int i, j;
	int neighbours=0;
	int local_rows, local_columns;

#ifdef TERMINATION_CHECK
	int	not_duplicate=1, not_dead=1;
	int	global_duplicate, global_dead;
#endif

	/*coords*/
	int coords[2];
	int neighbour_coords[2];

	/*neighbour ranks*/
	int north_rank, south_rank;

	/*requests for Isend/IRcv*/
	MPI_Request ISReqs[2], IRReqs[2];
	MPI_Status ISStatus[2], IRStatus[2];

	/*variables used for calculating time*/
	double local_start, local_finish, local_elapsed, elapsed;
	double local_overall_start, local_overall_finish, local_overall_elapsed, overall_elapsed;
	int n = 0;

	/*constructing mpi space*/
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*setting values for MPI_Cart_create*/
	ndims = 2;			/*  2D matrix/grneighbour_id */
	dim_size[0] = size;	/* rows */
	dim_size[1] = 1;	/* columns */
	periods[0] = 1;			/* row periodic */
	periods[1] = 1;			/* column periodic */
	reorder = 1;			/* allows processes reordered for efficiency */
	
	/*constructing cartesian matrix for processes(topology)*/
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &comm);

	/*we use cart_coords so we can know the coordinates of the process*/
	MPI_Cart_coords(comm, rank, ndims, coords);

	/*constructing the cell array & next phase cell array*/
	local_rows = ROW_SIZE/sqrt(size) + 2;
	local_columns = COLUMN_SIZE;
			
	/*make a new datatype so we can pass rows to other processes easier*/
	MPI_Datatype row;
	MPI_Type_vector(local_columns, 1, local_columns, MPI_INT, &row);
	MPI_Type_commit(&row);
	
	if( (cells = malloc(local_rows * local_columns * sizeof(int))) == NULL )
		perror_exit("Malloc failed(cells)");

	if( (np_cells = malloc( local_rows * local_columns * sizeof(int))) == NULL )
		perror_exit("Malloc failed(np_cells)");

	srand(time(NULL) + rank);

	for(i=1; i<local_rows-1; i++) {
		for(j=1; j<local_columns-1; j++) {
			cells[local_columns*i +j] = rand() % 2;
		}
	}

	/*using cart rank we found the neighbours of our process*/
	
	/*north neighbour*/
	neighbour_coords[0] = coords[0];
	neighbour_coords[1] = coords[1];
	neighbour_coords[0] -= 1;
	MPI_Cart_rank(comm, neighbour_coords, &north_rank);

	/*south neighbour*/
	neighbour_coords[0] = coords[0];
	neighbour_coords[1] = coords[1];
	neighbour_coords[0] += 1;
	MPI_Cart_rank(comm, neighbour_coords, &south_rank);

	/*synchronize processes*/
	MPI_Barrier(comm);
	/*================================================> start overall calculation time*/
	local_overall_start = MPI_Wtime();
	while( n < MAX_TIMES ) {
		/*================================================> start calculation time*/
		local_start = MPI_Wtime();

		/*send to the neighbours the appropriate columns and rows(non-blocking)*/
		//MPI_Isend( const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request)
		MPI_Isend(&cells[local_columns], local_columns, MPI_INT, north_rank, 0, comm, &ISReqs[0]);
		MPI_Isend(&cells[local_columns * (local_rows - 2)], local_columns, MPI_INT, south_rank, 1, comm, &ISReqs[1]);

		/*receive from the neighbours the appropriate columns and rows*/
		//int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
		MPI_Irecv(&cells[0], local_columns, MPI_INT, north_rank, 1, comm, &IRReqs[0]);
		MPI_Irecv(&cells[local_columns * (local_rows - 1)], local_columns, MPI_INT, south_rank, 0, comm, &IRReqs[1]);

		/*calculate the next phase without the info from the neighbours*/
		for(i=2; i<local_rows-2; i++) {
			for(j=0; j<local_columns; j++) {
				neighbours = Calculate_Neighbours(cells, local_columns, i, j);
				np_cells[local_columns * i + j] = Dead_Or_Alive(cells, local_columns, i, j, neighbours);

				neighbours=0;
			}
		}
		
		/*wait until we receive and send everything from/to neighbours*/
		MPI_Waitall(2, ISReqs, ISStatus);
		MPI_Waitall(2, IRReqs, IRStatus);

		/*continue calculating with the info from the neighbours we received*/

		/*calculating first and last row*/
		for(j=0; j<local_columns; j++) {
			neighbours = Calculate_Neighbours(cells, local_columns, 1, j);
			np_cells[local_columns + j] = Dead_Or_Alive(cells, local_columns, 1, j, neighbours);

			neighbours = 0;
		}

		for(j=0; j<local_columns; j++) {
			neighbours = Calculate_Neighbours(cells, local_columns, local_rows-2, j);
			np_cells[local_columns*(local_rows - 2) + j] = Dead_Or_Alive(cells, local_columns, local_rows-2, j, neighbours);

			neighbours = 0;
		}
		

		/*calculating first and last column*/
		for(i=1; i<local_rows-1; i++) {
			neighbours = Calculate_Neighbours(cells, local_columns, i, 0);
			np_cells[local_columns*i + 1] = Dead_Or_Alive(cells, local_columns, i, 1, neighbours);

			neighbours = 0;
		}

		for(i=1; i<local_rows-1; i++) {
			neighbours = Calculate_Neighbours(cells, local_columns, i, local_columns-1);
			np_cells[local_columns*i + local_columns - 1] = Dead_Or_Alive(cells, local_columns, i, local_columns - 1, neighbours);

			neighbours = 0;
		}

#ifdef TERMINATION_CHECK
		if( n % TERMCHECK_TIMES == 0) {	
			not_duplicate = 0;
			not_dead = 0;

			/*compare current phase with the next one && check if everything is dead*/
			for(i=1; i<local_rows-1; i++) {
				for(j=0; j<local_columns; j++) {
					if( cells[local_columns*i + j] != np_cells[local_columns*i + j] )
						not_duplicate = 1;
				}
			}
			
			for(i=1; i<local_rows-1; i++) {
				for(j=0; j<local_columns; j++) {
					if( cells[local_columns*i + j] == ALIVE )
						not_dead = 1;
				}
			}

			/*All reduce to check if everything is dead*/
			MPI_Allreduce(&not_dead, &global_dead, 1, MPI_INT, MPI_SUM, comm);
			if( global_dead == 0 ) {
				if( rank == 0)
					printf("Processes about to exit because every cell is dead!\n");

				free(cells);
				cells = NULL;

				free(np_cells);
				np_cells = NULL;

				MPI_Finalize();

				exit(EXIT_SUCCESS);
			}

			
			/*All reduce to check if next generation is the same as the current one*/
			MPI_Allreduce(&not_duplicate, &global_duplicate, 1, MPI_INT, MPI_SUM, comm);
			if( global_duplicate == 0 ) {
				if( rank == 0 )
					printf("Processes about to exit because the next generation is the same as the current one!\n");

				free(cells);
				cells = NULL;

				free(np_cells);
				np_cells = NULL;

				MPI_Finalize();

				exit(EXIT_SUCCESS);
			}
		}
#endif

		/*===================================================> end calculation time*/
		local_finish = MPI_Wtime();
		local_elapsed = local_finish - local_start;		

		if (n % REPEAT_TIMES == 0) {
			MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
			local_elapsed=0;

			if (rank == 0)
				printf("->Elapsed time = %.10f seconds\n", elapsed);
		}

		/*swap next generation array with the current generation array*/
		temp = cells;
		cells = np_cells;
		np_cells = temp;
		
		/*increment counter for loop */
		n++;
	} 
	/*================================================> end overall calculation time*/
	local_overall_finish = MPI_Wtime();
	local_overall_elapsed = local_overall_finish - local_overall_start;

	MPI_Reduce(&local_overall_elapsed, &overall_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	if(rank==0)
		printf("->>Total time elapsed = %.10f seconds\n" , overall_elapsed); 

	/*free allocated memory*/
	free(cells);
	cells = NULL;

	free(np_cells);
	np_cells = NULL;

	MPI_Finalize();

	exit(EXIT_SUCCESS);
}
