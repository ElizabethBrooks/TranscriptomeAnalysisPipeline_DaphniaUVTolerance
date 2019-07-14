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
repCount=0
treCount=0
genCount=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve inputs for number of reads, replicates, genotypes, and treatments
readMax=${#READARRAY[@]}-1
repMax=${#REPARRAY[@]}-1
treMax=${#TREARRAY[@]}-1
genMax=${#GENARRAY[@]}-1
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	#Determine what analysis method was used for the input folder of data
	if [[ $f1 == *"hisat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod=hisat2
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "$f1"/out/*.bam; do
			while [ $COUNTER -lt $readMax ]; do
				#Determine which read to add next to the set of replicates/samples
				if [[ $f2 == *"${REPARRAY[repCounter]}_${TREARRAY[genCounter]}_${GENARRAY[treCounter]}"* ]]; then
					if [[ $repCounter -eq $repMax && $treCounter -ne $treMax && $genCounter -ne $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2"
						let COUNTER+=1
						repCounter=0
						let treCounter+=1
					elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -ne $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2"
						let COUNTER+=1
						repCounter=0
						treCounter=0
						let genCounter+=1
					elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -eq $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2"
						let COUNTER+=1
					else
						#Add the next sample to the read array for input to cuffdiff
						READARRAY[COUNTER]="$f2, "
						let COUNTER+=1
						let repCounter+=1
					fi
				fi
			done
		done
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod=tophat2	
		#Loop through all forward and reverse paired reads and store the file locations in an array
		for f2 in "$f1"/out/*.bam; do
			while [ $COUNTER -lt $readMax ]; do
				#Determine which read to add next to the set of replicates/samples
				if [[ $f2 == *"${REPARRAY[repCounter]}_${TREARRAY[genCounter]}_${GENARRAY[treCounter]}"* ]]; then
					if [[ $repCounter -eq $repMax && $treCounter -ne $treMax && $genCounter -ne $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2/accepted_hits.bam"
						let COUNTER+=1
						repCounter=0
						let treCounter+=1
					elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -ne $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2/accepted_hits.bam"
						let COUNTER+=1
						repCounter=0
						treCounter=0
						let genCounter+=1
					elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -eq $genMax ]]; then
						#Add the last sample to the end of the set of replicates/samples
						READARRAY[COUNTER]="$f2/accepted_hits.bam"
						let COUNTER+=1
					else
						#Add the next sample to the read array for input to cuffdiff
						READARRAY[COUNTER]="$f2/accepted_hits.bam, "
						let COUNTER+=1
						let repCounter+=1
					fi
				fi
			done
		done
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
	echo ${READARRAY[@]}
	#Run cuffdiff on the aligned reads stored in the file array using 8 threads
	#cuffdiff -p 8 -o stats_"$analysisMethod"Tuxedo_run"$runNum" "$genomeFile" "${READARRAY[@]}"
done
