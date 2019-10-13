#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_cuffdiff_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/cufflinks/2.2.1
#Prepare for analysis
cd ..
dirFlag=0
runNum=0
COUNTER=0
repCount=0
treCount=0
genCount=0
readFlag=0
analysisTag=".sorted.bam"
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve input genome file path
inputGenomeFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePath.txt"
genomeFile=$(head -n 1 $inputGenomeFile)
#Retrieve inputs for number of reads, replicates, genotypes, and treatments
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/statsInputs_tuxedo.txt"
while IFS= read -r line; do
	#for word in $line; do
	#Each line contains the tags for the replicates, genotypes, or treatments
	#with each tag for the category separated by spaces
	for word in $line; do
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
		   	echo "ERROR: Incorrect number of lines in statsInputs_tuxedo.txt... exiting"
		   	exit 1
		fi
	done
	#done	
	let COUNTER+=1
done < "$inputsFile"
COUNTER=0
#Retrieve the number of replicates, genotypes, and samples
repMax=${#REPARRAY[@]}-1
treMax=${#TREARRAY[@]}-1
genMax=${#GENARRAY[@]}-1
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
		analysisMethod="hisat2"
	elif [[ $f1 == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"	
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Make a new directory for each analysis run
	while [ $dirFlag -eq 0 ]; do
		outputFolder=counts_cuffdiff_run"$runNum"
		mkdir "$outputFolder"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of Cuffdiff counting of $f1 data..."
		fi
	done
	#Loop through all forward and reverse paired reads and store the file locations in an array
	while [ $COUNTER -lt $readMax ]; do
		for f2 in "$f1"/*; do
			#Determine which read to add next to the set of replicates/samples
			if [[ $f2 == *${REPARRAY[repCounter]}"_"${GENARRAY[genCounter]}"_"${TREARRAY[treCounter]}* ]]; then
				if [[ $COUNTER -eq $readMax-1 ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY+="$f2$analysisExtension"
					LABELARRAY+="${GENARRAY[genCounter]}_${TREARRAY[treCounter]}"
					let COUNTER+=1					
				elif [[ $repCounter -eq $repMax && $treCounter -ne $treMax && $genCounter -ne $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY+="$f2$analysisExtension "
					LABELARRAY+="${GENARRAY[genCounter]}_${TREARRAY[treCounter]},"
					let COUNTER+=1
					repCounter=0
					let treCounter+=1
				elif [[ $repCounter -eq $repMax && $treCounter -ne $treMax && $genCounter -eq $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY+="$f2$analysisExtension "
					LABELARRAY+="${GENARRAY[genCounter]}_${TREARRAY[treCounter]},"
					let COUNTER+=1
					repCounter=0
					let treCounter+=1
				elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -ne $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY+="$f2$analysisExtension "
					LABELARRAY+="${GENARRAY[genCounter]}_${TREARRAY[treCounter]},"
					let COUNTER+=1
					repCounter=0
					treCounter=0
					let genCounter+=1
				elif [[ $repCounter -eq $repMax && $treCounter -eq $treMax && $genCounter -eq $genMax ]]; then
					#Add the last sample to the end of the set of replicates/samples
					READARRAY+="$f2$analysisExtension "
					LABELARRAY+="${GENARRAY[genCounter]}_${TREARRAY[treCounter]},"
					let COUNTER+=1
				else
					#Add the next sample to the read array for input to cuffdiff
					READARRAY+="$f2$analysisExtension,"
					let COUNTER+=1
					let repCounter+=1
				fi	
			fi
		done
	done
	#Double check that all input files were found
	#based on the number of reads specified in the inputsFile
	if [[ ${#READARRAY[@]} -ne $readMax ]]; then
		echo "ERROR: The number of reads identified for analysis does not match statsInputs_tuxedo... exiting"
		exit 1
	fi
	#Run cuffdiff on the aligned reads stored in the file array using 8 threads
	echo "Beginning statistical analysis of the following data set: "
	echo ${READARRAY[@]}
	echo "The following labels will be used to identify samples: "
	echo ${LABELARRAY[@]}
	cuffdiff -p 8 -L ${LABELARRAY[@]} -o "$outputFolder" "$genomeFile" ${READARRAY[@]}
	echo "Statistical analysis complete!"
done
