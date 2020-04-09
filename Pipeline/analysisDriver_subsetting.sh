#!/bin/bash
#Script to perform specified analysis of all samples in an input set
#Usage: bash analysisDriver.sh analysisMethod analysisArgs sampleList
#Usage Ex: bash analysisDriver.sh decoding Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of analysis method and sample names
if [ $# -eq 0 ]; then
   	echo "ERROR: No analysis method or sample name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=0
if [ "$1" == assemblyGenome ]; then #These are analysis methods that require additional input args
	#Loop through all input sets of treatments and perform t-test analsysis
	for i in "$@"; do
		#Skip first two arguments
		if [ $counter -ge 2 ]; then
			qsub assembly_genomeGuided_trinity.sh sortedCoordinate_samtoolsHisat2_run1 "$i" "$2"
		fi
		counter=$(($counter+1))
	done
else #These are analysis methods that do not require additional input args
	#Loop through all input sets of treatments and perform t-test analsysis
	for i in "$@"; do
		#Skip first argument
		if [ $counter -ge 1 ]; then
			if [ "$1" == decoding ]; then
				qsub decoding_transdecoder.sh trimmed_run1"$i"_assemblyTrinity
			elif [ "$1" == assembly ]; then
				qsub assembly_trinity.sh trimmed_run1 "$i"
			else
				echo "ERROR: Input analysis method is not a valid option... exiting!"
				exit 1
			fi
		fi
		counter=$(($counter+1))
	done
fi