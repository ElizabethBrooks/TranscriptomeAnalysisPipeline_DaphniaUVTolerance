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
runNum=1
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
		analysisMethod="Hisat2"
	elif [[ $f1 == *"tophat2"* ]]; then
		echo "ERROR: Tophat aligned files do not need to be sorted... exiting"
		exit 1
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		outputFolder=sorted_samtools"$analysisMethod"_run"$runNum"
		mkdir "$outputFolder"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for run $runNum of Samtools sorting of $f1 data..."
		fi
	done
	#Name output file of inputs
	inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
	#Sort input bam files if folder does not already exist
	if [ $? -eq 0 ]; then
		echo "Creating folder for sorted bam files..."
		#Loop through all reads and sort bam files for input to samtools
		for f2 in "$f1"/*/; do
			#Name of aligned file
			curAlignedSample="$f2"/accepted_hits.bam
			#Trim extension from current file name
			curSample=$(echo $f2 | sed 's/\.bam//')
			#Trim file path from current file name
			curSampleNoPath=$(basename $f2)
			curSampleNoPath=$(echo $curSampleNoPath | sed 's/\.bam//')
			#Create directory for current sample outputs
			mkdir "$outputFolder"/"$curSampleNoPath"
			#Run samtools to prepare mapped reads for sorting by name
			#using 8 threads
			echo "Sample $curSampleNoPath is being sorted..."
			samtools sort -@ 8 -n -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$curAlignedSample"
			echo "Sample $curSampleNoPath has been sorted!"
			#Add run inputs to output summary file
			echo "$curSampleNoPath" >> $inputOutFile
			echo samtools sort -@ 8 -n -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$curAlignedSample" >> $inputOutFile
		done
		#Copy previous summaries
		cp "$f1"/*.txt "$outputFolder"
	else
		echo "Sorted files already exists, skipping sorting..."
	fi
done