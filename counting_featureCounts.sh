#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_featureCounts_jobOutput
#$ -pe smp 8
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
analysisTag=".sorted.bam"
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve inputs for gff absolute path
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt"
genomeFile=$(head -n 1 $inputsFile)
#Retrieve outputs absolute path
outputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/outputsPath.txt"
outputsPath=$(head -n 1 $outputsFile)
#Move to outputs directory
cd "$outputsPath"
#Retrieve folders to analyze from the input arguments to the script
for f1 in $@; do
	#Determine if the folder name was input in the correct format
	if [[ $f1 == *\/* ]] || [[ $f1 == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi
	#Determine if the correct analysis folder was input
	if [[ $f1  != sorted* ]]; then
		echo "ERROR: The $f1 folder of aligned bam files were not found... exiting"
		exit 1
	fi
	#Determine what analysis method was used for the input folder of data
	if [[ $f1 == *"isat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="Hisat2"
		#Determine if hisat2 files were sorted
		if [[ $f1 == *"sorted"* ]]; then
			echo "Sorted Hisat2 files found! Proceeding with counting..."
		else
			echo "ERROR: The $f1 folder of bam files need to be sorted... exiting"
			exit 1
		fi
	elif [[ $f1 == *"ophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="Tophat2"
		#Determine if hisat2 files were sorted
		if [[ $f1 == *"sorted"* ]]; then
			echo "Sorted Tophat2 files found! Proceeding with counting..."
		else
			echo "ERROR: The $f1 folder of bam files need to be sorted... exiting"
			exit 1
		fi
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		outputFolder=counted_featureCounts"$analysisMethod"_run"$runNum"
		mkdir "$outputFolder"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for run $runNum of featureCounts counting of $f1 data..."
		fi
	done
	#Name output file of inputs
	inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
	#Loop through all sorted forward and reverse paired reads and store the file locations in an array
	for f2 in "$f1"/*/; do
		#Name of aligned file
		curAlignedSample="$f2"accepted_hits.bam
		#Trim file path from current file name
		curSampleNoPath=$(basename $f2)
		curSampleNoPath=$(echo $curSampleNoPath | sed 's/\.bam//')
		#Create directory for current sample outputs
		mkdir "$outputFolder"/"$curSampleNoPath"
		#Count reads using featureCounts
		echo "Sample $curSampleNoPath is being counted..."
		featureCounts -T 8 -p -t gene -g ID -a "$genomeFile" -s 2 --donotsort -o "$outputFolder"/"$curSampleNoPath"/counted.sam -M "$curAlignedSample"
		echo "Sample $curSampleNoPath has been counted!"
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> $inputOutFile
		echo featureCounts -T 8 -p -t gene -g ID -a "$genomeFile" -s 2 --donotsort -o "$outputFolder"/"$curSampleNoPath"/counted.sam -M "$curAlignedSample" >> $inputOutFile
	done
	#Copy previous summaries
	cp "$f1"/*.txt "$outputFolder"
done
