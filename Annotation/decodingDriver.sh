#!/bin/bash
#Script to perform specified analysis of all samples in an input set
#Usage: bash decodingDriver.sh analysisMethod analysisArgs sampleList
#Usage Ex: bash decodingDriver.sh decodingPB trimmed_run1 ncbi Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash decodingDriver.sh decodingPB sortedCoordinate_samtoolsHisat2_run1 ncbi Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash decodingDriver.sh decoding trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash decodingDriver.sh decoding sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of analysis method, data folder, and sample name(s)
if [ $# -lt 3 ]; then
   	echo "ERROR: No analysis method or sample name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=1
if [[ "$1" == decoding ]]; then #These are analysis methods that require one additional input args
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
			echo "Running decoding_transdecoder.sh $inputFolder"
			#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder
			qsub decoding_transdecoder.sh "$inputFolder"
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == decodingPB ]]; then #These are analysis methods that require two additional input args
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
			echo "Running decoding_transdecoderPfamBlastp.sh $inputFolder $3"
			#Usage: qsub decoding_transdecoderPfamBlastp.sh deNovoAssembledTranscriptomeFolder databaseSelection
			qsub decoding_transdecoderPfamBlastp.sh "$inputFolder" "$3"
		fi
		counter=$(($counter+1))
	done
else
	echo "ERROR: Input analysis method is not a valid option... exiting!"
	exit 1
fi