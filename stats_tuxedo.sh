#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_tuxedo_jobOutput
#$ -pe smp 8

#Prepare for alignment
cd ..
dirFlag=0
runNum=0
#Make a new directory for each alignment run
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
module load bio/cufflinks/2.2.1
COUNTER=0
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in aligned_tophat2/out/*; do
        READARRAY[COUNTER]="$f1/accepted_hits.bam, "
        let COUNTER+=1
done
#Re set the last array element to rmove the last two characters
# (extra comma and white space) from the last element of the read file array
unset 'READARRAY[${#READARRAY[@]}-1]'
READARRAY[COUNTER]="$f1/accepted_hits.bam"
genomeFile=TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff
#Run cuffdiff on the aligned reads stored in the file array using 8 threads
cuffdiff -p 8 -o stats_tuxedo_run"$runNum" "$genomeFile" "${READARRAY[@]}"