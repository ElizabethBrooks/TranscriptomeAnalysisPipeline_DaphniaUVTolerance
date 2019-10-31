#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sorting_samtools_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/python/2.7.14
module load bio/htseq/0.11.2
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
COUNTER=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve folders to analyze from the input arguments to the script
for f1 in "$@"; do
	#Determine if the folder name was input in the correct format
	if [[ $f1 == *\/* ]] || [[ $f1 == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi	
	#Determine what analysis method was used for the input folder of data
	if [[ $f1 == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
	elif [[ $f1 == *"tophat2"* ]]; then
		echo "ERROR: Tophat aligned files do not need to be sorted... exiting"
		exit 1
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		outputFolder=sorted_"$analysisMethod"Samtools_run"$runNum"
		mkdir "$outputFolder"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of Samtools sorting of $f1 data..."
		fi
	done
	#Sort input bam files if folder does not already exist
	if [ $? -eq 0 ]; then
		echo "Creating folder for sorted bam files..."
		#Loop through all reads and sort bam files for input to samtools
		for f3 in "$f1"/out/*; do
			#Trim extension from current file name
			curFile=$(echo $f3 | sed 's/\.bam//')
			#Trim file path from current file name
			curFileNoPath=$(basename $f3)
			curFileNoPath=$(echo $curFileNoPath | sed 's/\.bam//')
			echo "Sample $curFileNoPath is being sorted..."
			#Run samtools to prepare mapped reads for sorting by name
			#using 8 threads
			samtools sort -@ 8 -n -o "$outputFolder/$curFileNoPath".sorted.bam -T /tmp/"$analysisMethod"_sorted_"$f3".sorted "$f3"
			echo "Sample $curFileNoPath has been sorted!"
		done
	else
		echo "Sorted files already exists, skipping sorting..."
	fi
done