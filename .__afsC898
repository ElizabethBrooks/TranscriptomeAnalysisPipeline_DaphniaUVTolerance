#!/bin/bash
#Move out of script repository folder
cd ..
dirFlag=0
runNum=0
#Initializing file tag variables
genotypeTag=""
treatmentTag=""
replicationTag=""
#Initialize counters and flags
wordCOUNTER=0
fileFlag=0
#Retrieve location of gene count files to be merged
geneCounts=$1
#Determine if the folder name was input in the correct format
#if [[ $geneCounts == *\/* ]] || [[ $geneCounts == *\\* ]]; then
	#echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	#exit 1
#fi	
#Set name for merge gene counts file
mergedCounts="geneCounts_merged.txt"
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	outputFolder=counts_merged_run"$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of merging counts of $1 data..."
	fi
done
#Retrieve merge order list from file
#UPDATE: similarly to statsInputs_tuxedo.txt for counting_cuffdiff.sh
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/mergeInputs.txt"
#Merge gene counts from file based on order in mergeOrder.txt
while IFS= read -r line; do
	#Determine tags for current file in list
	for word in $line; do
		if [[ wordCOUNTER -eq 0 ]]; then
			REPARRAY[repCount]="$word"
		   	let repCount+=1
		elif [[ wordCOUNTER -eq 1 ]]; then
		   	TREARRAY[treCount]="$word"
		   	let treCount+=1
		elif [[ wordCOUNTER -eq 2 ]]; then
		   	GENARRAY[genCount]="$word"
		   	let genCount+=1
		fi
	done
	let wordCOUNTER+=1
done < "$inputsFile"
wordCOUNTER=0
#Create array of sorted tags for search files
for genTag in ${GENARRAY[@]}; do
	for treTag in ${TREARRAY[@]}; do
		for repTag in ${REPARRAY[@]}; do
			TAGARRAY[wordCOUNTER]="$repTag"_"$genTag"_"$treTag"
			let wordCOUNTER+=1
		done
	done
done
wordCOUNTER=0
#Merge files based on tag order
for currentFile in ${TAGARRAY[@]}; do
	echo "Sample $currentFile is being merged..."
	if [ $fileFlag -eq 0 ]; then #Output the first column with gene IDs
		cp "$geneCounts"/*"$currentFile"* "$outputFolder"/"$mergedCounts"
		#Insert header line
		sed -i.bak 1i"gene0" "$outputFolder"/"$mergedCounts"
		fileFlag=1
	else #Add the gene counts from the next file
		echo $(paste -d' ' "$outputFolder"/"$mergedCounts" <(cut -d' ' -f1 "$geneCounts"/*"$currentFile"*)) >> "$outputFolder"/"$mergedCounts"
	fi
	#Insert current file tags to header line
	sed -i.bak "1 s/$/ $currentFile/" "$outputFolder"/"$mergedCounts"
	let wordCOUNTER=0
done