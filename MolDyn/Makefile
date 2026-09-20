CC     = gcc
CFLAGS = -O2 -fopenmp
CHUNK ?= 0

all: md_serial md_omp

md_serial: main.c md_util.c forces.c moldyn.h
	$(CC) $(CFLAGS) main.c md_util.c forces.c -o md_serial

md_omp: main.c md_util.c forces_omp.c moldyn.h
	$(CC) $(CFLAGS) -DCHUNK=$(CHUNK) main.c md_util.c forces_omp.c -o md_omp

clean:
	rm -f md_serial md_omp *.o
