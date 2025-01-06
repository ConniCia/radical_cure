#!/bin/bash

# Command line input
PATH_ON_CLUSTER=$1
CURRENT_ARRAY_FILE=$2
STEP=$3


# Get parameters for current tasks from array file passed as argument
echo "Working inside job array file ${CURRENT_ARRAY_FILE} at line $((SLURM_ARRAY_TASK_ID+1))"
echo ""


# Read parameter names
NAMES=$(sed "1q;d" ${CURRENT_ARRAY_FILE})
NAMES_VECTOR=($NAMES)


# Read parameter values of CURRENT job
VALUES=$(sed "$((SLURM_ARRAY_TASK_ID+1))q;d" ${CURRENT_ARRAY_FILE})
VALUES_VECTOR=($VALUES)


# Assign job parameters
PARS_GEN=$PATH_ON_CLUSTER${VALUES_VECTOR[0]}
PARS_DOM=$PATH_ON_CLUSTER${VALUES_VECTOR[1]}
PARS_OCC=$PATH_ON_CLUSTER${VALUES_VECTOR[1]} # same file for domestic and occupational mosquitoes
PARS_INT=$PATH_ON_CLUSTER${VALUES_VECTOR[2]}
PARS_OUT=$PATH_ON_CLUSTER${VALUES_VECTOR[3]}
SAMPLE=${VALUES_VECTOR[4]}
SEED=1 # seed has to be at least 1 in Michael's code


# create output directory if it doesn't yet exist
DIR_OUT=$(dirname ${PARS_OUT})
mkdir -p $DIR_OUT


# print looped parameters to file
if [ "$SAMPLE" -eq 1 ];
then
  echo ${NAMES_VECTOR[@]} > ${DIR_OUT}/row_varying.tsv
  echo ${VALUES_VECTOR[@]} >> ${DIR_OUT}/row_varying.tsv
fi


# Print task info to task output file
echo "Calling Pv model with the following arguments:"
echo -e "\tSeed: ${SEED}"
echo -e "\tGeneral parameters: ${PARS_GEN}"
echo -e "\tDomestic mosquitos parameters: ${PARS_DOM}"
echo -e "\tOccupational mosquitos parameters: ${PARS_OCC}"
echo -e "\tInterventions: ${PARS_INT}"
echo -e "\tOutput: ${PARS_OUT}"
echo ""
echo -e "Calling: ${PATH_ON_CLUSTER}/data/${STEP}/input/RP13_executable ${SEED} ${PARS_GEN} ${PARS_DOM} ${PARS_OCC} ${PARS_INT} ${PARS_OUT}"


# srun the task
srun --job-name=individual --partition=common,dedicated --qos=fast --mem-per-cpu=300 ${PATH_ON_CLUSTER}/data/${STEP}/input/RP13_executable ${SEED} ${PARS_GEN} ${PARS_DOM} ${PARS_OCC} ${PARS_INT} ${PARS_OUT}
