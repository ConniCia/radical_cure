#!/bin/bash
module load R/4.2.1

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2
N_LHS=$3
N_MCMC=$4
N_CHAINS=$5
BURNIN=$6

# Select current input file
INPUTS=$(ls ${PATH_ON_CLUSTER}/data/${STEP}/input/input*.csv)
VECTOR_INPUTS=($INPUTS)
CURRENT_FILE=${VECTOR_INPUTS[${SLURM_ARRAY_TASK_ID} - 1]}

# Run Rscript
Rscript ${PATH_ON_CLUSTER}/code/R/_run_fit.R ${PATH_ON_CLUSTER} ${STEP} ${CURRENT_FILE} ${N_LHS} ${N_MCMC} ${N_CHAINS} ${BURNIN}