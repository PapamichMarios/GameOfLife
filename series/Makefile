CC = mpicc -O3
EXE = mpi.x
FIL = mpi.c
CFLAGS = -Wall -g
MATH = -lm

all: $(EXE)

mpi.x: mpi.o
	$(CC) $(CFLAGS) -o mpi.x mpi.c ../funcs.c $(MATH)

clean:
	rm -f $(EXE)

count:
	wc $(FIL)

