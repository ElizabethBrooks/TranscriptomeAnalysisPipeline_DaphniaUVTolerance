#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash exactTestDriver_edgeR.sh sampleList
#Usage Ex: bash exactTestDriver_edgeR.sh Y05 Y023_5 E05 R2 PA Sierra

#Load module for R
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Skip first argument
	if [ $counter -ge 1 ]; then
		bash exactTest_subset_edgeR.sh $i
	fi
	counter=$(($counter+1))
done