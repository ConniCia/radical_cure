#!/bin/bash
module load R/4.2.1

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2


# List output directories
OUTPUT_DIRS=$(ls -d ${PATH_ON_CLUSTER}/data/${STEP}/output/*/)
VECTOR_OUTPUT_DIRS=($OUTPUT_DIRS)
CURRENT_DIR=${VECTOR_OUTPUT_DIRS[${SLURM_ARRAY_TASK_ID} - 1]}

echo $CURRENT_DIR


# Run Rscript
Rscript ${PATH_ON_CLUSTER}/code/R/_run_clean.R ${PATH_ON_CLUSTER} ${CURRENT_DIR}
