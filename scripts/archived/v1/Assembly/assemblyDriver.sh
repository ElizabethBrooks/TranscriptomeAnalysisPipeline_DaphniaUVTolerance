#!/bin/bash
#Script to perform specified analysis of all samples in an input set
#Usage: bash assemblyDriver.sh analysisMethod analysisArgs sampleList
#Usage Ex: bash assemblyDriver.sh genomeAssembly sortedCoordinate_samtoolsHisat2_run1 14239 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash assemblyDriver.sh assembly trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of analysis method, data folder, and sample name(s)
if [ $# -lt 3 ]; then
   	echo "ERROR: No analysis method or sample name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=1
if [[ "$1" == assembly ]]; then #These are analysis methods that require one additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Determine what type of data folder was input
		if [[ "$2" == trimmed* ]]; then
			inputFolder=$(echo "$2""$i"_assemblyTrinity)
		elif [[ "$2" == sorted* ]]; then
			inputFolder=$(echo "$2""$i"_assemblyGenomeTrinity)
		else
			echo "ERROR: Input folder for analysis is not a valid option... exiting!"
			exit 1
		fi
		#Skip first two arguments
		if [ $counter -ge 3 ]; then
			echo "Running assembly_trinity.sh $2 $i"
			#Usage: qsub assembly_trinity.sh trimmedFolder genotype
			qsub assembly_trinity.sh "$2" "$i"
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == genomeAssembly ]]; then #These are analysis methods that require two additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Determine what type of data folder was input
		if [[ "$2" == trimmed* ]]; then
			inputFolder=$(echo "$2""$i"_assemblyTrinity)
		elif [[ "$2" == sorted* ]]; then
			inputFolder=$(echo "$2""$i"_assemblyGenomeTrinity)
		else
			echo "ERROR: Input folder for analysis is not a valid option... exiting!"
			exit 1
		fi
		#Skip first three arguments
		if [ $counter -ge 4 ]; then
			echo "Running assembly_genomeGuided_trinity.sh $2 $i $3"
			#Usage: qsub assembly_genomeGuided_trinity.sh sortedFolder genotype maxIntronLength
			qsub assembly_genomeGuided_trinity.sh "$2" "$i" "$3"
		fi
		counter=$(($counter+1))
	done
else
	echo "ERROR: Input analysis method is not a valid option... exiting!"
	exit 1
fi