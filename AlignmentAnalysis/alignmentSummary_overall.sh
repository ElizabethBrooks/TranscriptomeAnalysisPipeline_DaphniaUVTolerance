#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_overall.sh alignmentFolder
#Usage Ex: bash alignmentSummary_overall.sh aligned_tophat2_run2
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve folders to analyze from the input arguments to the script
for f1 in $@; do
	#Determine if the folder name was input in the correct format
	if [[ "$f1" == *\/* ]] || [[ "$f1" == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi
	#Determine what analysis method was used for the input folder of data
	if [[ "$f1" == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
	elif [[ "$f1" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Prepare input and output file names
	inputStats=../../"$f1"
	outputStats=../../alignmentSummarized_"$analysisMethod"
	#Retrieve run number for input alignment folder
	runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
	#Retrieve summaries for each aligned sample
	for f2 in "$inputStats"/*/; do
		echo "Merging sample $f2 of $analysisMethod alignment summary..."
		#Retrieve sample summary based on alignment method
		bash alignmentSummary_"$analysisMethod"_sample.sh "$f1" "$analysisMethod"
		#Combine summaries into one csv file
		cat "$outputStats"_combined_"$runNum".csv >> "$outputStats"_allSamples_"$runNum".csv
		rm "$outputStats"_combined_"$runNum".csv
		echo "Sample $f2 of $analysisMethod alignment summary has been merged!"
	done
done