#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "timer.h"
#include "funcs.h"

#define ROW_SIZE	3840		
#define	COLUMN_SIZE	3840

#define REPEAT_TIMES 	20
#define MAX_TIMES	500
#define TERMCHECK_TIMES 10

// #define TERMINATION_CHECK

int main() {

	int * cells = NULL; 
	int * np_cells = NULL;
	int * temp = NULL;
	int	i, j;
	int local_rows, local_columns;
	double start, loop_start;
	double finish, loop_finish;	
	int n = 0;
	int y, Up, Down, Left, Right, neighbours=0;

#ifdef TERMINATION_CHECK
	int not_dead=0;
	int not_duplicate=0;
#endif

	local_columns = COLUMN_SIZE;
	local_rows = ROW_SIZE;
	
	if( (cells = malloc(local_rows * local_columns * sizeof(int))) == NULL )
	       perror_exit("Malloc(1) failed");
		
	if( (np_cells = malloc( local_rows * local_columns * sizeof(int))) == NULL )
	       perror_exit("Malloc(2) failed");

	srand(time(NULL));
	for(i=0; i<local_rows; i++) {
		for(j=0; j<local_columns; j++) {
			cells[local_columns*i +j] = rand() % 2;
			np_cells[local_columns*i +j] = 0;
		}
	}

#ifdef DEBUG
	for(i=0; i<ROW_SIZE; i++) {
		for(j=0; j<COLUMN_SIZE; j++) {
			printf("%d ", cells[COLUMN_SIZE*i +j]);
		}
		printf("\n");
	}
#endif

	/*================================================> start overall calculation time*/
	GET_TIME(start);
	while( n < MAX_TIMES ) {

		/*================================================> start loop calculation time*/
		GET_TIME(loop_start);
	  	for (i = 0; i < local_rows; ++i) {

			Up = ((i + local_rows - 1) % local_rows) * local_columns;
			y = i * local_columns;
			Down = ((i + 1) % local_rows) * local_columns;

			for ( j = 0; j < local_columns; ++j) {
		 		Left = (j + local_columns- 1) % local_columns;
		 		Right = (j + 1) % local_columns;

				neighbours =  cells[Left + Up] 
					    + cells[j + Up]
					    + cells[Right + Up]
					    + cells[Left + y]
					    + cells[Right + y]
					    + cells[Left + Down] 
					    + cells[j + Down] 
					    + cells[Right + Down];
					
				if( (neighbours == 3) || ( (cells[i*local_columns + j] ==  ALIVE) && (neighbours == 2) ) )
					np_cells[i*local_columns+j] = ALIVE;
				else
					np_cells[i*local_columns+j] = DEAD;
			}
		}

#ifdef TERMINATION_CHECK
		if( n % TERMCHECK_TIMES == 0) {	
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

				exit(EXIT_SUCCESS);
			}
		}
#endif

		/*================================================> end loop calculation time*/
		GET_TIME(loop_finish);
		if( n % REPEAT_TIMES == 0 )
			printf("->Elapsed loop time = %.10f seconds\n", loop_finish-loop_start);

		/*swap next generation array with the current generation array*/
		temp = cells;
		cells = np_cells;
		np_cells = temp;
		
		/*increment counter for loop */
		n++;
	}
	/*================================================> end overall calculation time*/
	GET_TIME(finish);

#ifdef DEBUG
	printf("\n");
	for(i=0; i<ROW_SIZE; i++) {
		for(j=0; j<COLUMN_SIZE; j++) {
			printf("%d ", cells[COLUMN_SIZE*i +j]);
		}
		printf("\n");
	}
#endif
	printf("->>Elapsed overall time = %.10f seconds\n" ,  finish-start);

	/*free allocated memory*/
	free(cells);
	cells = NULL;

	free(np_cells);
	np_cells = NULL;

	exit(EXIT_SUCCESS);
}
