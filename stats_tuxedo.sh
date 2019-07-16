#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_tuxedo_jobOutput
#$ -pe smp 8 
#Required modules for ND CRC servers
module load bio/cufflinks/2.2.1
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
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/statsInputs_tuxedo.txt"
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
	    	echo "Incorrect number of lines in statsInputs_tuxedo... exiting"
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
		analysisMethod=hisat2
		analysisTag=
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod=tophat2	
		analysisTag="/accepted_hits.bam"
	else
		echo "The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Loop through all forward and reverse paired reads and store the file locations in an array
	while [ $COUNTER -lt $readMax ]; do
		for f2 in "$f1"/out/*; do
			#Determine which read to add next to the set of replicates/samples
			if [[ $f2 == *${REPARRAY[repCounter]}"_"${GENARRAY[genCounter]}"_"${TREARRAY[treCounter]}* ]]; then
				if [[ $COUNTER -eq $readMax-1 ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY[COUNTER]="$f2$analysisTag"
					let COUNTER+=1					
				elif [[ $repCounter -eq $repMax && $treCounter -ne $treMax && $genCounter -ne $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY[COUNTER]="$f2$analysisTag"
					let COUNTER+=1
					repCounter=0
					let treCounter+=1
				elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -ne $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY[COUNTER]="$f2$analysisTag"
					let COUNTER+=1
					repCounter=0
					treCounter=0
					let genCounter+=1
				elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -eq $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY[COUNTER]="$f2$analysisTag"
					let COUNTER+=1
				else
					#Add the next sample to the read array for input to cuffdiff
					READARRAY[COUNTER]="$f2$analysisTag, "
					let COUNTER+=1
					let repCounter+=1
				fi	
			fi
		done
	done
	#Double check that all input files were found
	#based on the number of reads specified in the inputsFile
	if [[ ${#READARRAY[@]} -ne $readMax ]]; then
		echo "The number of reads identified for analysis does not match statsInputs_tuxedo... exiting"
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
	#echo ${READARRAY[@]}
	#Loop through all reads and sort bam files for input to cuffdiff
	for f2 in ${#READARRAY[@]}; do
		echo "Sample ${f2:13:${#f2}-3} is being sorted..."
		#Run samtools to prepare mapped reads for sorting
		#using 8 threads
		samtools sort -@ 8 -o stats_tuxedo_run"$runNum"/${f2:13:${#f2}-3}.sorted.bam -T /tmp/${f2:13:${#f2}-3}.sorted $f2
		echo "Sample ${f2:13:${#f2}-3} has been sorted!"
	done
	#Run cuffdiff on the aligned reads stored in the file array using 8 threads
	cuffdiff -p 8 -o stats_"$analysisMethod"Tuxedo_run"$runNum" "$genomeFile" ${READARRAY[@]}
done
