#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput
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
analysisTag=".sorted.bam"
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve inputs for gff absolute path
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt"
genomeFile=$(head -n 1 $inputsFile)
#Retrieve folders to analyze from the input arguments to the script
for f1 in "$@"; do
	#Determine if the folder name was input in the correct format
	if [[ $f1 == *\/* ]] || [[ $f1 == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		outputFolder=counted_htseq_run"$runNum"
		mkdir "$outputFolder"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for run $runNum of HtSeq counting of $f1 data..."
		fi
	done
	#Name output file of inputs
	inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
	#Loop through all sorted forward and reverse paired reads and store the file locations in an array
	for f2 in "$f1"/*; do
		#Name of aligned file
		curAlignedSample="$f2"/accepted_hits.bam
		#Trim file path from current file name
		curSampleNoPath=$(basename $f2)
		curSampleNoPath=$(echo $curSampleNoPath | sed 's/\.bam//')
		echo "Sample $curSampleNoPath is being counted..."
		#Create directory for current sample outputs
		mkdir "$outputFolder"/"$curSampleNoPath"
		#Sort sam file reads by name for counting
		samtools sort -n -T /tmp/"$outputFolder"/"$curSampleNoPath"/sorted.bam -o "$outputFolder"/"$curSampleNoPath"/sorted.bam "$curAlignedSample"
		#Count reads using htseq-count
		htseq-count -f bam -s no -m union -t gene -i ID -o "$outputFolder"/"$curSampleNoPath"/counted.sam "$outputFolder"/"$curSampleNoPath"/sorted.bam "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		echo "Sample $curSampleNoPath has been counted!"
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> $inputOutFile
		echo samtools sort -n -T /tmp/"$outputFolder"/"$curSampleNoPath"/sorted.bam -o "$outputFolder"/"$curSampleNoPath"/sorted.bam "$curAlignedSample" >> $inputOutFile
		echo htseq-count -f bam -s no -m union -t gene -i ID -o "$outputFolder"/"$curSampleNoPath"/counted.sam "$outputFolder"/"$curSampleNoPath"/sorted.bam "$genomeFile" ">" "$outputFolder"/"$curSampleNoPath"/counts.txt >> $inputOutFile
	done
done
