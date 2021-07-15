#!/usr/bin/env bash
CURR_DIR=`pwd`
BASE="$CURR_DIR/.."
HEAD="-I$BASE/inc"
gcc -c -O3 $HEAD $BASE/runs/gen_data.c
gcc -c -O3 $HEAD $BASE/src/multilayerD.c
gcc -c -O3 $HEAD $BASE/src/activation.c
gcc -c -O3 $HEAD $BASE/src/ran2.c

gcc -o $BASE/bin/gen_data.out gen_data.o activation.o multilayerD.o  ran2.o -O3 -lm 

rm *.o
