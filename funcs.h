#ifndef __FUNCS_H__
#define __FUNCS_H__

#define	ALIVE 1
#define	DEAD  0


void perror_exit(char* message);
int Calculate_Neighbours(const int * cells, const int local_columns, const int i, const int j);
int Dead_Or_Alive(const int * cells, const int local_columns, const int i, const int j, const int neighbours);

#endif