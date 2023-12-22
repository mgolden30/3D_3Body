CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -llapacke -O3


main: main.c clsqr2/lsqr.o
	$(CC) $^ -o $@ $(CFLAGS)

#compilation command found in local README
clsqr2/lsqr.o: clsqr2/lsqr.c
	cd clsqr2; $(CC)  -pedantic -Wall -c -I. lsqr.c; cd ..

clean:
	rm main clsqr2/lsqr.o


