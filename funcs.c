#include <stdlib.h>
#include <stdio.h>

#include "funcs.h"

void perror_exit(char * message) {
	perror(message);
	exit(EXIT_FAILURE);
}


int Calculate_Neighbours(const int * cells, const int local_columns, const int i, const int j) {
	int neighbours = cells[local_columns * (i-1) + j] + 	/*north neighbour*/
			 cells[local_columns * (i+1) + j] + 	/*south neighbour*/
			 cells[local_columns * i     + (j+1)] + /*east neighbour*/
			 cells[local_columns * i     + (j-1)] + /*west neighbour*/
			 cells[local_columns * (i-1) + (j-1)] + /*north west neighbour*/
			 cells[local_columns * (i-1) + (j+1)] + /*north east neighbour*/
			 cells[local_columns * (i+1) + (j-1)] + /*south west cell neighbour*/
			 cells[local_columns * (i+1) + (j+1)];  /*south east cell neighbour*/ 

	return neighbours;
}

int Dead_Or_Alive(const int * cells, const int local_columns, const int i, const int j, const int neighbours) {
	/*-if an array cell has 3 neighbours, it is born or it remains alive
	  -if an array cell is alive and it has 2 exact neighbours it remains alive
	
	  else it dies or remains dead
	*/
	
	if( (neighbours == 3) || ( (cells[local_columns * i +j] == ALIVE) && (neighbours == 2) ) )
		return ALIVE;

	return DEAD;
}