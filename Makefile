CC=gcc
CFLAGS=-lm -lgsl -lgslcblas -O3


main: main.c clsqr2/lsqr.o
	$(CC) $^ -o $@ $(CFLAGS)

#compilation command found in local README
clsqr2/lsqr.o: clsqr2/lsqr.c
	cd clsqr2; $(CC)  -pedantic -Wall -c -I. lsqr.c; cd ..

#Automate downloading the LSQR files
clsqr2/lsqr.c:
	wget https://web.stanford.edu/group/SOL/software/lsqr/c/clsqr2.zip; \
	unzip clsqr2.zip; \
	rm clsqr2.zip;


clean:
	rm main; \
	rm clsqr2 -r


