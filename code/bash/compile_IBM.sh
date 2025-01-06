#!/bin/bash
module load gcc/9.2.0

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2


# compile
g++ -O3 -o ${PATH_ON_CLUSTER}/data/${STEP}/input/RP13_executable ${PATH_ON_CLUSTER}/code/Pv_mod/Pv_mod/*.cpp
