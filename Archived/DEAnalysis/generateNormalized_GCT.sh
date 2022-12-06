#!/bin/bash
#Script to perform reformatting of normalized count sets to GCT format
#Usage: bash reformatCountsDriver_normalized_GCT.sh countsFolder sampleList
#Usage Ex: bash reformatCountsDriver_normalized_GCT.sh GeneCountsAnalyzed_countedCoordinate_htseqHisat2_run1_fullset Y05 Y023_5 E05 R2 PA Sierra
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
		bash reformatCounts_normalized_GCT.sh "$1" $i
	fi
	counter=$(($counter+1))
done