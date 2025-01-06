#!/bin/bash

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2


# Run array of jobs: run a job for each distinct parameter combination
ARRAY_FILE=(${PATH_ON_CLUSTER}/data/${STEP}/input/varying.tsv)
N_LINES=$(wc -l < ${ARRAY_FILE})

if [ $N_LINES -gt 1 ]; # only run jobs if there are jobs to run... duh!
then
  JOB_ID=$(sbatch --job-name=RP13_${STEP} --partition=common,dedicated --qos=fast --mem-per-cpu=300 --array=1-${N_LINES} ${PATH_ON_CLUSTER}/code/bash/run_IBM_individual.sh ${PATH_ON_CLUSTER} ${ARRAY_FILE} ${STEP} | sed 's/Submitted batch job //')
  echo -e "Submitted batch job ${JOB_ID}: running ${N_LINES} simulations"
fi
