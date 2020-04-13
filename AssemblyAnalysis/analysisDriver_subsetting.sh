#!/bin/bash
#Script to perform specified alignment of all assembly folders in an input set
#Usage: bash analysisDriver_subsetting.sh alignmentSoftware assemblyFolderList
#Usage Ex: bash analysisDriver_subsetting.sh hisat2 trimmed_run1E05_assemblyTrinity trimmed_run1PA_assemblyTrinity trimmed_run1R2_assemblyTrinity trimmed_run1Sierra_assemblyTrinity trimmed_run1Y023_5_assemblyTrinity trimmed_run1Y05_assemblyTrinity

#Check for input arguments of analysis method and sample names
if [ $# -eq 0 ]; then
   	echo "ERROR: No alignment method or assembly folder name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=0
#Loop through all input assembly folders
for i in "$@"; do
	#Skip first argument
	if [ $counter -ge 1 ]; then
		if [ "$1" == "hisat2" ]; then
			qsub ../Pipeline/alignment_hisat2.sh "$i"
		elif [ "$1" == "tophat2" ]; then
			qsub ../Pipeline/alignment_tophat2.sh "$i"
		else
			echo "ERROR: Input alignment method is not a valid option... exiting!"
			exit 1
		fi
	fi
	counter=$(($counter+1))
done