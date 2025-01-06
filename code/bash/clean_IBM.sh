#!/bin/bash

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2


# Count the number of output directories
N_DIRS=$(ls -d ${PATH_ON_CLUSTER}/data/${STEP}/output/*/ | wc -l)


# Run R script to clean output files in each output directory
JOB_ID1=$(sbatch --job-name=RP13_clean --partition=common,dedicated --qos=fast --array=1-${N_DIRS} ${PATH_ON_CLUSTER}/code/bash/clean_IBM_individual.sh ${PATH_ON_CLUSTER} ${STEP} | sed 's/Submitted batch job //')
echo -e "Submitted batch job ${JOB_ID1}: cleaning IBM output from ${N_DIRS} folders"


# Run bash script to concatenate csv files created in above sbatch
JOB_ID2=$(sbatch --job-name=RP13_concatenate --dependency=afterok:$JOB_ID1 --partition=common,dedicated --qos=fast ${PATH_ON_CLUSTER}/code/bash/clean_IBM_concatenate.sh ${PATH_ON_CLUSTER} ${STEP} | sed 's/Submitted batch job //')
echo -e "Submitted batch job ${JOB_ID2}: concatenating resulting csv files"
