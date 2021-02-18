#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash exactTestDriver_edgeR.sh sampleList
#Usage Ex: bash exactTestDriver_edgeR.sh Y05 Y023_5 E05 R2 PA Sierra

#Load module for R
#module load bio

#Check for input arguments of sample names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	bash exactTest_subset_edgeR.sh $i
done