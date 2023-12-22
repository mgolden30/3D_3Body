This is a periodic orbit finder for the 3 body problem in 3D. It uses two 
linear algbera libraries that can instaled in various ways.

1. GSL - the GNU Scientific Library for its implementation of cblas.
   If you have a better cblas implementation, by all means use it.

2. clsqr2 - the most recent implementation of the lsqr algorithm in C.
            Find the source code here: 
            C files 2: Contributed Aug 2007 by Michael Friedlander, UBC
            https://web.stanford.edu/group/SOL/software/lsqr/

            NOTE: The Makefile will compile this code for you! 
            Just put clsqr2/ in this directory

To compile
	make main

To uninstall 
	make clean



Just run ./main after compilation to start hunting orbits! States will be saved to the states/ directory in the format [r1; r2; p1; p2; T] where r's and p's are 3D vectors. Happy hunting!
