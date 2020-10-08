CC = mpicc
GCC = gcc
CE = mpiexec
BIN = mpi.o open_mp.o
EXE = mpi open_mp
FIL = mpi.c open_mp.c
CFLAGS = -Wall -g
MATH = -lm
OPENMP = -fopenmp

all: $(EXE)

mpi: mpi.o
	$(CC) $(CFLAGS) -o mpi mpi.o ../funcs.o $(MATH)

mpi.o: mpi.c
	$(CC) $(CFLAGS) -c mpi.c $(MATH)

open_mp: open_mp.o
	$(CC) $(CFLAGS) -o open_mp open_mp.o ../funcs.o $(MATH) $(OPENMP)

open_mp.o: open_mp.c
	$(CC) $(CFLAGS) -c open_mp.c $(MATH) $(OPENMP)

clean:
	rm -f $(BIN) $(EXE)

count:
	wc $(FIL)

