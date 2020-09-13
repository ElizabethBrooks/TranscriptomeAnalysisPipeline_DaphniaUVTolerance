#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash annotationDriver_edgeR.sh analysisType countsFolder sampleList
#Usage Ex: bash annotationDriver_edgeR.sh exactTest genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash annotationDriver_edgeR.sh 2WayANOVA_filtered genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1 Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Initialize variables
counter=0
#Loop through all input sets of treatments and perform t-test analsysis
for i in "$@"; do
	#Skip first two arguments
	if [ $counter -ge 2 ]; then
		bash geneAnnotations_edgeR.sh "$1" "$2" $i
	fi
	counter=$(($counter+1))
done