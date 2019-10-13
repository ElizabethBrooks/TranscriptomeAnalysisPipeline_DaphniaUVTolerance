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
#Retrieve merge order for file
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/mergeOrder.txt"
while IFS= read -r line; do
	for word in $line; do
		if [[ COUNTER -eq 0 ]]; then
			replicationTag=$word
		elif [[ COUNTER -eq 1 ]]; then
			genotypeTag="$word"
		elif [[ COUNTER -eq 2 ]]; then
		   	treatmentTag="$word"
		else
		   	echo "ERROR: Incorrect number of tags in line for mergeOrder.txt... exiting"
		   	exit 1
		fi
		let COUNTER+=1
		#Merge files based on tag order
		if [ $fileFlag -eq 0 ]; then #Output the first column with gene IDs
			cp "$geneCounts"/*"$replicationTag"_"$genotypeTag"_"$treatmentTag"* geneCounts_merged.txt
		else
			paste -d" " geneCounts_merged.txt "$geneCounts"/*"$replicationTag"_"$genotypeTag"_"$treatmentTag"*
		else
			echo "ERROR: Please check that mereOrder.txt input tags are in the same order found in gene count file names... exiting!"
		fi
	done
done < "$inputsFile"