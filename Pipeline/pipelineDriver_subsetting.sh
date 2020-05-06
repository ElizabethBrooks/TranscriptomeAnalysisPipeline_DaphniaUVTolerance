#!/bin/bash
#Script to perform specified analysis of all samples in an input set
#Usage: bash pipelineDriver_subsetting.sh analysisMethod analysisArgs sampleList
#Usage Ex: bash pipelineDriver_subsetting.sh decoding trimmed_run1 ncbi Y05 Y023_5 E05 R2 PA Sierra
#Alternate usage Ex: bash pipelineDriver_subsetting.sh genomeAssembly sortedCoordinate_samtoolsHisat2_run1 14239 Y05 Y023_5 E05 R2 PA Sierra
#Alternate usage Ex: bash pipelineDriver_subsetting.sh buildingHisat2 trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Alternate usage Ex: bash pipelineDriver_subsetting.sh buildingHisat2 sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Alternate usage Ex: bash pipelineDriver_subsetting.sh alignmentHisat2 trimmed_run1 20 14239 Y05 Y023_5 E05 R2 PA Sierra
#Alternate usage Ex: bash pipelineDriver_subsetting.sh alignmentHisat2 sortedCoordinate_samtoolsHisat2_run1 20 14239 Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of analysis method, data folder, and sample name(s)
if [ $# -lt 3 ]; then
   	echo "ERROR: No analysis method or sample name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=1
if [[ "$1" == assembly || "$1" == buildingTophat2 || "$1" == buildingHisat2 ]]; then #These are analysis methods that require one additional input args
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
			if [[ "$1" == assembly ]]; then
				#Usage: qsub assembly_trinity.sh trimmedFolder genotype
				qsub assembly_trinity.sh "$2" "$i"
			elif [[ "$1" == buildingTophat2 || "$1" == buildingBowtie2 ]]; then
				#Usage: qsub building_bowtie2.sh trimmedOrAssemblyFolder
				qsub building_bowtie2.sh "$inputFolder"
			elif [[ "$1" == buildingHisat2 ]]; then
				#Usage: qsub building_hisat2.sh trimmedOrAssemblyFolder
				bash building_hisat2.sh "$inputFolder"
			fi
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == genomeAssembly || "$1" == decoding ]]; then #These are analysis methods that require two additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Skip first three arguments
		if [ $counter -ge 4 ]; then
			if [ "$1" == genomeAssembly ]; then
				#Usage: qsub assembly_genomeGuided_trinity.sh sortedFolder genotype maxIntronLength
				qsub assembly_genomeGuided_trinity.sh "$2" "$i" "$3"
			elif [[ "$1" == decoding ]]; then
				#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder databaseSelection
				qsub decoding_transdecoder.sh "$inputFolder" "$3"
			fi
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == alignmentTophat2 || "$1" == alignmentHisat2 ]]; then #These are analysis methods that require three additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Skip first 4 arguments
		if [ $counter -ge 5 ]; then
			#Run slected alignment software
			if [[ "$1" == alignmentTophat2 ]]; then
				#Usage: qsub alignment_tophat2.sh trimmedOrAssemblyFolder minIntronLength maxIntronLength
				qsub alignment_tophat2.sh "$inputFolder" "$3" "$4"
			elif [[ "$1" == alignmentHisat2 ]]; then
				#Usage: qsub alignment_hisat2.sh trimmedOrAssemblyFolder minIntronLength maxIntronLength
				qsub alignment_hisat2.sh "$inputFolder" "$3" "$4"
			fi
		fi
		counter=$(($counter+1))
	done
else
	echo "ERROR: Input analysis method is not a valid option... exiting!"
	exit 1
fi