#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio/python/2.7.14
module load bio/htseq/0.11.2
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
genomeFile=TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	#Determine what analysis method was used for the input folder of data
	if [[ $f1 == *"hisat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod=hisat2
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "$f1"/out/*.bam; do
	    	READARRAY[COUNTER]="$f2, "
			let COUNTER+=1
		done
		#Re set the last array element to remove the last two characters
		unset 'READARRAY[${#READARRAY[@]}-1]'
		READARRAY[COUNTER]="$f2"
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod=tophat2	
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "$f1"/out/*; do
	    	READARRAY[COUNTER]="$f2/accepted_hits.bam, "
			let COUNTER+=1
		done
		#Re set the last array element to remove the last two characters
		unset 'READARRAY[${#READARRAY[@]}-1]'
		READARRAY[COUNTER]="$f2/accepted_hits.bam"
	else
		echo "The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		mkdir stats_"$analysisMethod"Tuxedo_run"$runNum"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of tuxedo stats analysis of $f1 data..."
			#Reset the folder name flag for different analysis methods
			let runNum=0
		fi
	done
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		mkdir stats_edgeR_run"$runNum"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of edgeR stats analysis..."
		fi
	done
	#Loop through all forward and reverse paired reads and store the file locations in arrays
	for f2 in "$f1"/*pForward.fq.gz; do
		echo "Sample ${f2:13:${#f2}-28} is being sorted and counted..."
		#Run samtools to prepare mapped reads for counting
		# using 8 threads
		samtools sort -@ 8 -o stats_edgeR_run"$runNum"/"${f2:13:${#f2}-28}"/accepted_hits.sorted.bam -T /tmp/"${f2:13:${#f2}-28}"/accepted_hits.sorted.bam $f1/accepted_hits.bam
		#Run htseq-count to prepare sorted reads for stats analysis in edgeR
		htseq-count -s no -m union -t gene -i trID $f1/accepted_hits.sorted.bam -i "$genomeFile" > stats_edgeR_run"$runNum"/"${f2:13:${#f2}-28}".counts
		echo "Sample ${f2:13:${#f2}-28} has been counted!"
	done
done
