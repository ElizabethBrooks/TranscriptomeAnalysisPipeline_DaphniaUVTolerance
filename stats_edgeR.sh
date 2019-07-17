#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#$ -pe smp 4
#Required modules for ND CRC servers
module load bio
module load bio/python/2.7.14
module load bio/htseq/0.11.2
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
COUNTER=0
genomeFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff"
repCount=0
treCount=0
genCount=0
readFlag=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve inputs for number of reads, replicates, genotypes, and treatments
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/statsInputs_edgeR.txt"
while IFS= read -r line; do
	for word in $line; do
		#Each line contains the tags for the replicates, genotypes, or treatments
		#with each tag for the category separated by spaces
	    if [[ COUNTER -eq 0 ]]; then
	    	readMax=$word
	    elif [[ COUNTER -eq 1 ]]; then
	    	REPARRAY[repCount]="$word"
	    	let repCount+=1
	    elif [[ COUNTER -eq 2 ]]; then
	    	TREARRAY[treCount]="$word"
	    	let treCount+=1
	    elif [[ COUNTER -eq 3 ]]; then
	    	GENARRAY[genCount]="$word"
	    	let genCount+=1
	    else
	    	echo "Incorrect number of lines in statsInputs_edgeR... exiting"
	    	#exit 1
	    fi
	done	
	let COUNTER+=1
done < "$inputsFile"
#Retrieve the number of replicates, genotypes, and samples
repMax=${#REPARRAY[@]}-1
treMax=${#TREARRAY[@]}-1
genMax=${#GENARRAY[@]}-1
#Retrieve folders to analyze from the input arguments to the script
COUNTER=0
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
		analysisTag=""
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"	
		analysisTag="/accepted_hits.bam"
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
			#Reset the folder name flag for different analysis methods
			let runNum=0
		fi
	done
	#IF
	#Loop through all reads and sort bam files for input to cuffdiff
	for f3 in "$f1"/out/*; do
		echo "Sample ${f3:24:${#f3}-(28+${#analysisTag})} is being sorted..."
		#Run samtools to prepare mapped reads for sorting
		#using 4 threads
		samtools sort -@ 4 -o stats_"$analysisMethod"EdgeR_run"$runNum"/${f3:24:${#f3}-(28+${#analysisTag})}.sorted.bam -T /tmp/"$analysisMethod"EdgeR/${f3:24:${#f3}-(28+${#analysisTag})}.sorted $f3
		echo "Sample ${f3:24:${#f3}-(28+${#analysisTag})} has been sorted!"
	done
	#Loop through all forward and reverse paired reads and store the file locations in an array
	for f2 in stats_"$analysisMethod"EdgeR_run"$runNum"/*.sorted.bam; do
		READARRAY[COUNTER]="$f2, "
		let COUNTER+=1					
	done
	unset 'READARRAY[COUNTER-1]'
	READARRAY[COUNTER-1]="$f2"
	#Run htseq-count to prepare sorted reads for stats analysis in edgeR
	echo "Beginning statistical analysis of the following data set: "
	echo ${READARRAY[@]}
	#ORDER
	htseq-count -s no -m union -t gene -i trID -o stats_"$analysisMethod"EdgeR_run"$runNum"/out.counted.sam ${READARRAY[@]} -i "$genomeFile"
	echo "Reads have been counted!"
done
