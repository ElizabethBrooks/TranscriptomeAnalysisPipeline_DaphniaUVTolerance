#!/bin/bash
#Initializing file tag variables
genotypeTag=""
treatmentTag=""
replicationTag=""
#Initialize counters and flags
COUNTER=0
fileFlag=0
#Retrieve location of gene count files to be merged
geneCounts=$(head -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/geneCountsPath.txt")
#Set name for merge gene counts file
mergedCounts="geneCounts_merged.txt"
#Retrieve merge order list from file
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/mergeOrder.txt"
#Merge gene counts from file based on order in mergeOrder.txt
while IFS= read -r line; do
	#Determine tags for current file in list
	for word in $line; do
		if [[ COUNTER -eq 0 ]]; then
			replicationTag=$word
		elif [[ COUNTER -eq 1 ]]; then
			genotypeTag="$word"
		elif [[ COUNTER -eq 2 ]]; then
		   	treatmentTag="$word"
		else
		   	echo "ERROR: Incorrect number of tags in a line for mergeOrder.txt... exiting"
		   	exit 1
		fi
		let COUNTER+=1
	done
	#Merge files based on tag order
	currentFile="$replicationTag"_"$genotypeTag"_"$treatmentTag"
	if [ $fileFlag -eq 0 ]; then #Output the first column with gene IDs
		cp "$geneCounts"/*"$currentFile"* "$mergedCounts"
	else #Add the gene counts from the next file
		cut -d' ' -f1 "$geneCounts"/*"$currentFile"*
		paste -d' ' "$mergedCounts" "$geneCounts"/*"$currentFile"*
	else
		echo "ERROR: Please check that mereOrder.txt input tags are in the same order found in gene count file names... exiting!"
	fi
	#Add current file tags to header line
	sed -i.bak "1 s/$/ $currentFile/" "$mergedCounts"
	let COUNTER=0
done < "$inputsFile"