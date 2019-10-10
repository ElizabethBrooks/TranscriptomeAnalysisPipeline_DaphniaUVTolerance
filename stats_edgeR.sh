#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
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
#Retrieve inputs for gff absolute path
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePath.txt"
genomeFile=$(head -n 1 $inputsFile)
#Retrieve folders to analyze from the input arguments to the script
for f1 in "$@"; do
	#Determine if the folder name was input in the correct format
	if [[ $f1 == *\/* ]] || [[ $f1 == *\\* ]]; then
		echo "Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi	
	#Determine what analysis method was used for the input folder of data
	if [[ $f1 == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
		analysisTag=".bam"
		analysisFiles="stats_"$analysisMethod"EdgeR_sorted"
		analysisExtension=""
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"	
		analysisTag="/accepted_hits.bam"
		analysisFiles="$f1/out/"
		analysisExtension="/accepted_hits.bam"
	else
		echo "The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		mkdir stats_"$analysisMethod"EdgeR_run"$runNum"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of edgeR stats analysis of $f1 data..."
		fi
	done
	#Sort input bam files if folder does not already exist
	if [ "$analysisMethod" == "hisat2" ]; then
		mkdir "$analysisFiles"
		if [ $? -eq 0 ]; then
			echo "Creating folder for sorted bam files..."
			#Loop through all reads and sort bam files for input to cuffdiff
			for f3 in "$f1"/out/*; do
				echo "Sample ${f3:(${#f1}+5):(${#f3}-${#analysisTag}-${#analysisFiles})} is being sorted..."
				#Run samtools to prepare mapped reads for sorting by name
				#using 8 threads
				samtools sort -@ 8 -n -o "$analysisFiles/${f3:(${#f1}+5):(${#f3}-${#analysisTag}-${#analysisFiles})}".sorted.bam -T /tmp/"$analysisMethod"EdgeR_sorted_"${f3:(${#f1}+5):(${#f3}-${#analysisTag}-${#analysisFiles})}".sorted $f3
				echo "Sample ${f3:(${#f1}+5):(${#f3}-${#analysisTag}-${#analysisFiles})} has been sorted!"
			done
		else
			echo "Sorted files already exists, skipping sorting..."
		fi
	fi
	#Loop through all forward and reverse paired reads and store the file locations in an array
	for f2 in "$analysisFiles"/*; do
		echo "Sample $f2$analysisExtension is being counted..."
		htseq-count -f bam -s no -m union -t gene -i ID -o "stats_"$analysisMethod"EdgeR_run"$runNum"/${f2:${#analysisFiles}:(${#f2}-${#analysisFiles}+${#analysisExtension}+${#analysisTag}-15)}.counted.sam" "$f2$analysisExtension" "$genomeFile" > "stats_"$analysisMethod"EdgeR_run"$runNum"/${f2:${#analysisFiles}:(${#f2}-${#analysisFiles}+${#analysisExtension}+${#analysisTag}-15)}_counts.txt"
		echo "Sample $f2$analysisExtension has been counted!"
	done
done
