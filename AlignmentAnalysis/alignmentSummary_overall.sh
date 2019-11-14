#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_overall.sh alignmentFolder
#Usage Ex: bash alignmentSummary_overall.sh aligned_topaht2_run2
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Move to directory with output alignment folders
cd ../..
#Retrieve folders to analyze from the input arguments to the script
for f1 in $@; do
	folderInput=$(echo "$f1")
	#Determine if the folder name was input in the correct format
	if [[ "$folderInput" == *\/* ]] || [[ "$folderInput" == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi
	#Determine what analysis method was used for the input folder of data
	if [[ "$folderInput" == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
	elif [[ "$folderInput" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
	else
		echo "ERROR: The $folderInput folder or bam files were not found... exiting"
		#exit 1
	fi
	#Prepare input and output file names
	inputStats="$folderInput"
	outputStats=TranscriptomeAnalysisPipeline_DaphniaUVTolerance/AlignmentAnalysis/alignmentSummarized_"$analysisMethod"
	#Retrieve run number for input alignment folder
	runNum=$(echo "$folderInput" | sed "s/aligned_"$analysisMethod"_//g")
	#Retrieve summaries for each aligned sample
	for f2 in "$folderInput"/*/; do
		sampleInput=$(echo "$f2")
		echo "Merging sample $sampleInput of $analysisMethod alignment summary..."
		#Retrieve sample summary based on alignment method
		#bash alignmentSummary_"$analysisMethod"_sample.sh "$f1" "$analysisMethod"
		#Combine summaries into one csv file
		#cat "$outputStats"_"$analysisMethod"_combined_"$runNum".csv >> "$outputStats"_"$analysisMethod"_allSamples_"$runNum".csv
		echo "Sample $sampleInput of $analysisMethod alignment summary has been merged!"
	done
done