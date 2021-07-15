#!/usr/bin/env bash
CURR_DIR=`pwd`
BASE="$CURR_DIR/.."
HEAD="-I$BASE/inc"
gcc -c -O3 $HEAD $BASE/runs/generic_pot.c
gcc -c -O3 $HEAD $BASE/src/multilayerD.c
gcc -c -O3 $HEAD $BASE/src/simAnnealing.c
gcc -c -O3 $HEAD $BASE/src/adam.c
gcc -c -O3 $HEAD $BASE/src/minimize.c
gcc -c -O3 $HEAD $BASE/src/activation.c
gcc -c -O3 $HEAD $BASE/src/ran2.c

gcc -o $BASE/bin/generic_pot.out generic_pot.o simAnnealing.o activation.o minimize.o multilayerD.o adam.o ran2.o -O3 -lm -lgsl -lgslcblas

rm *.o
