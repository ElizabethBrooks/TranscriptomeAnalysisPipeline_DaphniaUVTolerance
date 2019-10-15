#!/bin/bash
#Move out of script repository folder
cd ..
#Initializing file tag variables
genotypeTag=""
treatmentTag=""
replicationTag=""
#Initialize counters and flags
wordCOUNTER=0
fileFlag=0
#Retrieve location of gene count files to be merged
geneCounts=$1
#Set name for merge gene counts file
mergedCounts="geneCounts_merged.txt"
mkdir "$mergedCounts"
#Retrieve merge order list from file
#UPDATE: similarly to statsInputs_tuxedo.txt for counting_cuffdiff.sh
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/mergeInputs.txt"
#Merge gene counts from file based on order in mergeOrder.txt
while IFS= read -r line; do
	#Determine tags for current file in list
	#UPDATE: add line tags to identify fields
	for word in $line; do
		if [[ COUNTER -eq 0 ]]; then
			replicationTag=$word
		elif [[ COUNTER -eq 1 ]]; then
			genotypeTag="$word"
		elif [[ COUNTER -eq 2 ]]; then
		   	treatmentTag="$word"
		else
		   	echo "ERROR: Incorrect number of lines in mergeInputs.txt... exiting"
		   	exit 1
		fi
		let wordCOUNTER+=1
	done
	#Merge files based on tag order
	echo "File $currentFile is being merged..."
	currentFile="$replicationTag"_"$genotypeTag"_"$treatmentTag"
	if [ $fileFlag -eq 0 ]; then #Output the first column with gene IDs
		cp "$geneCounts"/*"$currentFile"* "$mergedCounts"
		#Insert header line
		sed -i.bak 1i"gene0" "$mergedCounts"
	else #Add the gene counts from the next file
		cut -d' ' -f1 "$geneCounts"/*"$currentFile"*
		paste -d' ' "$mergedCounts" "$geneCounts"/*"$currentFile"*
	#else
		#echo "ERROR: Please check that mereOrder.txt input tags are in the same order found in gene count file names... exiting!"
	fi
	#Insert current file tags to header line
	sed -i.bak "1 s/$/ $currentFile/" "$mergedCounts"
	let wordCOUNTER=0
done < "$inputsFile"