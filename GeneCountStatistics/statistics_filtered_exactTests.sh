#!/bin/bash
#Script to perform t-test analysis of all samples in an input set
#Usage: bash statistics_exactTests.sh countsFile.csv maxCols
#Usage Ex: bash statistics_exactTests.sh ../GeneCounts_Merged/geneCounts_merged_counted_htseqTophat2_run1_fullset_cleaned.csv 36
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Prepare for analysis
iMax=0
colMax=$(($2-5))
#Loop through all sets of treatments and perform t-test analsysis
for i in $(seq 1 $colMax); do
	#Offset by sample treatment set
	if [[ $i -gt $iMax ]]; then
		iMax=$(($i+5))
		bash statistics_filtered_DEAnalysis_edgeR.sh "$1" $i $iMax
		i=$iMax
	fi
done