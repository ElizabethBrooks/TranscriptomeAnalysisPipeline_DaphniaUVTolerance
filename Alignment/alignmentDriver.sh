#!/bin/bash
#Script to perform specified analysis of all samples in an input set
#Usage: bash alignmentDriver.sh analysisMethod analysisArgs sampleList
#Usage Ex: bash alignmentDriver.sh buildingHisat2 trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash alignmentDriver.sh buildingHisat2 sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash alignmentDriver.sh alignmentHisat2 trimmed_run1 trimmed_run1 20 14239 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash alignmentDriver.sh alignmentHisat2 trimmed_run1 sortedCoordinate_samtoolsHisat2_run1 20 14239 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash alignmentDriver.sh sorting name aligned_hisat2_run2 sortedCoordinate_samtoolsTophat2_run1 Y05 Y023_5 E05 R2 PA Sierra
#Usage Ex: bash alignmentDriver.sh sorting name aligned_hisat2_run2 trimmed_run1 Y05 Y023_5 E05 R2 PA Sierra

#Check for input arguments of analysis method, data folder, and sample name(s)
if [ $# -lt 3 ]; then
   	echo "ERROR: No analysis method or sample name(s) supplied... exiting!"
   	exit 1
fi
#Initialize variables
counter=1
if [[ "$1" == buildingTophat2 || "$1" == buildingHisat2 ]]; then #These are analysis methods that require one additional input args
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
			if [[ "$1" == buildingTophat2 || "$1" == buildingBowtie2 ]]; then
				echo "Running building_bowtie2.sh $inputFolder"
				#Usage: qsub building_bowtie2.sh trimmedOrAssemblyFolder
				qsub building_bowtie2.sh "$inputFolder"
			elif [[ "$1" == buildingHisat2 ]]; then
				echo "Running building_hisat2.sh $inputFolder"
				#Usage: qsub building_hisat2.sh trimmedOrAssemblyFolder
				bash building_hisat2.sh "$inputFolder"
			fi
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == sorting ]]; then #These are analysis methods that require three additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Determine what type of data folder was input
		if [[ "$4" == trimmed* ]]; then
			inputFolder=$(echo "$4""$i"_assemblyTrinity)
		elif [[ "$4" == sorted* ]]; then
			inputFolder=$(echo "$4""$i"_assemblyGenomeTrinity)
		else
			echo "ERROR: Input folder for analysis is not a valid option... exiting!"
			exit 1
		fi
		#Skip first 4 arguments
		if [ $counter -ge 5 ]; then
			#Usage: qsub sorting_samtools.sh sortingTarget sortingMethod alignedFolder assembledFolder
			qsub sorting_samtools.sh assembly "$2" "$3" "$inputFolder"
		fi
		counter=$(($counter+1))
	done
elif [[ "$1" == alignmentTophat2 || "$1" == alignmentHisat2 ]]; then #These are analysis methods that require three additional input args
	#Loop through all input sets of treatments and perform selected analsysis
	for i in "$@"; do
		#Determine what type of data folder was input
		if [[ "$3" == trimmed* ]]; then
			inputFolder=$(echo "$3""$i"_assemblyTrinity)
		elif [[ "$3" == sorted* ]]; then
			inputFolder=$(echo "$3""$i"_assemblyGenomeTrinity)
		else
			echo "ERROR: Input folder for analysis is not a valid option... exiting!"
			exit 1
		fi
		#Skip first 5 arguments
		if [ $counter -ge 6 ]; then
			#Run slected alignment software
			if [[ "$1" == alignmentTophat2 ]]; then
				echo "Running alignment_tophat2.sh assembly $2 $inputFolder $4 $5"
				#Usage: qsub alignment_tophat2.sh alignmentTarget trimmedFolder assemblyFolder minIntronLength maxIntronLength
				qsub alignment_tophat2.sh assembly "$2" "$inputFolder" "$4" "$5"
			elif [[ "$1" == alignmentHisat2 ]]; then
				echo "Running alignment_hisat2.sh assembly $2 $inputFolder $4 $5"
				#Usage: qsub alignment_hisat2.sh alignmentTarget trimmedFolder assemblyFolder minIntronLength maxIntronLength
				qsub alignment_hisat2.sh assembly "$2" "$inputFolder" "$4" "$5"
			fi
		fi
		counter=$(($counter+1))
	done
else
	echo "ERROR: Input analysis method is not a valid option... exiting!"
	exit 1
fi