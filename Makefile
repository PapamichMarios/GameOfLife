GCC = gcc -O3
BIN = sequential.o
EXE = sequential
FIL = funcs.c sequential.c
CFLAGS = -Wall -g

.PHONY: all
all:
	make install
	cd blocks && make all
	cd series && make all

install: $(EXE)

sequential: 
	$(GCC) $(CFLAGS) -o sequential sequential.c funcs.c -lm

clean:
	rm -f $(BIN) $(EXE)
	cd blocks && make clean
	cd series && make clean

count:
	wc $(FIL)
