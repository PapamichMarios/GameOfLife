CC = mpicc -O3
EXE = mpi.x open_mp.x cuda.x
FIL = mpi.c open_mp.c cuda.cu
CFLAGS = -Wall -g
MATH = -lm
OPENMP = -fopenmp
CUDA = nvcc

all: $(EXE)

mpi.x:
	$(CC) $(CFLAGS) -o mpi.x mpi.c ../funcs.c $(MATH)

open_mp.x:
	$(CC) $(CFLAGS) -o open_mp.x open_mp.c ../funcs.c $(MATH) $(OPENMP)

cuda.x:
	$(CUDA) cuda.cu -o cuda.x ../funcs.c $(MATH)

clean:
	rm -f $(EXE)

count:
	wc $(FIL)

