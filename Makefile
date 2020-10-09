GCC = gcc
BIN = funcs.o sequential.o
EXE = sequential
FIL = funcs.c sequential.c
CFLAGS = -Wall -g

.PHONY: all
all:
	make install
	cd blocks && make all
	cd series && make all

install: funcs.o $(EXE)

sequential: sequential.o
	$(GCC) $(CFLAGS) -o sequential sequential.o funcs.o -lm

sequential.o: sequential.c
	$(GCC) $(CFLAGS) -c sequential.c

funcs.o: funcs.c
	$(GCC) $(CFLAGS) -c funcs.c -lm

clean:
	rm -f $(BIN) $(EXE)
	cd blocks && make clean
	cd series && make clean

count:
	wc $(FIL)
