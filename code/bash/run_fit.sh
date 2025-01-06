#!/bin/bash

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2
N_LHS=$3
N_MCMC=$4
N_CHAINS=$5
BURNIN=$6


# Each input file is a different job
N_INPUTS=$(ls ${PATH_ON_CLUSTER}/data/${STEP}/input/input*.csv | wc -l)


# Run R script to fit to IMPROV trial data
sbatch --job-name=RP13_${STEP} --partition=common,dedicated --qos=fast --array=1-${N_INPUTS} --cpus-per-task=${N_CHAINS} ${PATH_ON_CLUSTER}/code/bash/run_fit_individual.sh ${PATH_ON_CLUSTER} ${STEP} ${N_LHS} ${N_MCMC} ${N_CHAINS} ${BURNIN}
