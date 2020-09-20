#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "../timer.h"
#include "../funcs.h"

#define ROW_SIZE	360
#define	COLUMN_SIZE	360

#define REPEAT_TIMES 	20
#define MAX_TIMES	500	
#define TERMCHECK_TIMES 10

#define THREAD_SIZE	512

//#define TERMINATION_CHECK

__global__ void Game_of_Life_Kernel(int* cells, int * np_cells, int * size, int * columns) {
	int cell_index = threadIdx.x + blockIdx.x * blockDim.x;

	if( cell_index < *size )
	{ 
		/*find neighbours*/
		int x = cell_index % *columns;
		int y = cell_index - x;
		int Left  = (x + *(columns) - 1) % *columns;
		int Right = (x + 1) % *columns;
		int Up = (y + *size - *columns) % *size;
		int Down = (y + *columns) % *size;
		 
		int neighbours =  cells[Left + Up]     /*north west*/
				+ cells[x + Up]        /*north*/
				+ cells[Right + Up]    /*north east*/
				+ cells[Left + y]      /*west*/
				+ cells[Right + y]     /*east*/
				+ cells[Left + Down]   /*south west*/
				+ cells[x + Down]      /*south*/
				+ cells[Right + Down]; /*south east*/

		if( (neighbours == 3) || ( (cells[cell_index] ==  ALIVE) && (neighbours == 2) ) )
			np_cells[cell_index] = ALIVE;
		else
			np_cells[cell_index] = DEAD;
	}
}

/*---------------->*/
int main() {

	int * cells 	 = NULL;
	int * np_cells = NULL;
	int * temp = NULL;
	int array_size;
	int	array_columns;
	int	size;
	int i,j;
	int n = 0;
	double loop_start, loop_finish;
	double start, finish;

#ifdef TERMINATION_CHECK
	int	not_dead=0;
	int	not_duplicate=0;
#endif
	int * dcells = NULL;
	int * dnp_cells = NULL;
	int * darray_size = NULL;
	int * darray_columns = NULL;

	srand(time(NULL));
	size = ROW_SIZE * COLUMN_SIZE * sizeof(int);

	/*malloc the cell arrays*/
	if( (cells = (int*) malloc(size)) == NULL )
		exit(EXIT_FAILURE);

	if( (np_cells = (int*) malloc(size)) == NULL )
		exit(EXIT_FAILURE);

	/*how many threads in a block*/
	array_size = ROW_SIZE * COLUMN_SIZE;
	array_columns = COLUMN_SIZE;

	/*allocate device space*/
	cudaMalloc((void**)&dcells, size);
	cudaMalloc((void**)&dnp_cells, size);

	cudaMalloc((void**)&darray_size, sizeof(int));
	cudaMalloc((void**)&darray_columns, sizeof(int));

	/*setup input values*/
	for(i=0; i<ROW_SIZE; i++)
	{
		for(j=0; j<COLUMN_SIZE; j++)
		{
			cells[COLUMN_SIZE*i +j] = rand() % 2;
			np_cells[COLUMN_SIZE*i + j] = 0;
		}
	}

	/*================================================> start overall calculation time*/
	GET_TIME(start);
	while( n < MAX_TIMES)
	{
		/*================================================> start loop calculation time*/
		GET_TIME(loop_start);

		/*copy inputs to device*/
		cudaMemcpy(dcells, cells, size, cudaMemcpyHostToDevice);
		cudaMemcpy(dnp_cells, np_cells, size, cudaMemcpyHostToDevice);

		cudaMemcpy(darray_size, &array_size, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(darray_columns, &array_columns, sizeof(int), cudaMemcpyHostToDevice);

		/*call kernel*/
		Game_of_Life_Kernel<<<array_size/THREAD_SIZE, THREAD_SIZE>>>(dcells, dnp_cells, darray_size, darray_columns);

		/*copy result to host*/
		cudaMemcpy(cells, dcells, size, cudaMemcpyDeviceToHost);
		cudaMemcpy(np_cells, dnp_cells, size, cudaMemcpyDeviceToHost);

#ifdef TERMINATION_CHECK
		if( n % TERMCHECK_TIMES == 0)	
		{	
			not_duplicate = 0;
			not_dead = 0;

			/*compare current phase with the next one && check if everything is dead*/
			for(i=0; i<ROW_SIZE; i++) {
				for(j=0; j<COLUMN_SIZE; j++) {
					if( cells[COLUMN_SIZE*i + j] != np_cells[COLUMN_SIZE*i + j] )
						not_duplicate = 1;
				}
			}
			
			for(i=0; i<ROW_SIZE; i++) {
				for(j=0; j<COLUMN_SIZE; j++) {
					if( cells[COLUMN_SIZE*i + j] == ALIVE )
						not_dead = 1;
				}
			}

			if( not_dead == 0 )
				printf("!--->Every cell is dead, program about to exit.\n");

			if( not_duplicate == 0 )
				printf("!--->Current cell generation is the same as the next one.\n");

			if( not_dead == 0  || not_duplicate == 0) {
				/*free memory*/
				free(cells);
				free(np_cells);

				cells = NULL;
				np_cells = NULL;

				cudaFree(dcells);
				cudaFree(dnp_cells);
				cudaFree(darray_size);
				cudaFree(darray_columns);

				dcells = NULL;
				dnp_cells = NULL;
				darray_size = NULL;
				darray_columns = NULL;

				exit(EXIT_SUCCESS);
			}
		}
#endif
		/*================================================> finish loop calculation time*/
		GET_TIME(loop_finish);
		if(n % REPEAT_TIMES == 0)
			printf("->Elapsed time = %.10f seconds\n", loop_finish-loop_start);

		/*swap arrays*/
		temp = cells;
		cells = np_cells;
		np_cells = temp;

		/*increment loop counter*/
		n++;
	}

	/*================================================> finish overall calculation time*/
	GET_TIME(finish);
	printf("->>Elapsed overall time = %.10f seconds\n", finish-start);

	/*free memory*/
	free(cells);
	free(np_cells);

	cells = NULL;
	np_cells = NULL;

	cudaFree(dcells);
	cudaFree(dnp_cells);
	cudaFree(darray_size);
	cudaFree(darray_columns);

	dcells = NULL;
	dnp_cells = NULL;
	darray_size = NULL;
	darray_columns = NULL;

	exit(EXIT_SUCCESS);
}
