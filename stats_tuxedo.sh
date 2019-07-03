#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_tuxedo_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
#module load bio/cufflinks/2.2.1
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
COUNTER=0
genomeFile=TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	mkdir stats_tuxedo_run"$runNum"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of tuxedo stats analysis..."
	fi
done
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	if [[ $f1 == *"hisat2"* ]]; then
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "$f1"/out/*.bam; do
	    	READARRAY[COUNTER]="$f2.bam, "
			let COUNTER+=1
		done
		#Re set the last array element to rmove the last two characters
		# (extra comma and white space) from the last element of the read file array
		unset 'READARRAY[${#READARRAY[@]}-1]'
		READARRAY[COUNTER]="$f2.bam"
	elif [[ $f1 == *"tophat2"* ]]; then	
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "f1"/out/*; do
	    	READARRAY[COUNTER]="$f2/*.bam, "
			let COUNTER+=1
		done
		#Re set the last array element to rmove the last two characters
		# (extra comma and white space) from the last element of the read file array
		unset 'READARRAY[${#READARRAY[@]}-1]'
		READARRAY[COUNTER]="$f2/accepted_hits.bam"
	else
		echo "The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	echo ${READARRAY[@]}
	#Run cuffdiff on the aligned reads stored in the file array using 8 threads
	#cuffdiff -p 8 -o stats_tuxedo_run"$runNum" "$genomeFile" "${READARRAY[@]}"
done