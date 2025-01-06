#!/bin/bash

# Command line input
PATH_ON_CLUSTER=$1
STEP=$2


# create results folder if it doesn't yet exist
mkdir -p ${PATH_ON_CLUSTER}/data/${STEP}/results/


# Concatenate all cleaned csv files, but only print the header once!
for METRIC in {timeline,stats,effsize}; do
  FILE_RES="${PATH_ON_CLUSTER}/data/${STEP}/results/${METRIC}.csv"
  i=0                                            # reset counter
  for FILE_CLEAN in ${PATH_ON_CLUSTER}/data/${STEP}/results_partial/${METRIC}_*.csv; do
    echo $FILE_CLEAN
    if [[ $i -eq 0 ]] ; then
      head -1 "$FILE_CLEAN" > "$FILE_RES"      # copy header if it is the first file
    fi
      tail -n +2  "$FILE_CLEAN" >> "$FILE_RES" # for each other file, append from the 2nd line
    i=$(( $i + 1 ))                            # increment counter
    echo $i
  done
done
