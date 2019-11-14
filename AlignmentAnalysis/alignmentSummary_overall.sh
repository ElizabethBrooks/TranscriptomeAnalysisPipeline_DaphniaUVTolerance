#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_overall.sh alignmentFolder
#Usage Ex: bash alignmentSummary_overall.sh aligned_tophat2_run2
#Move to directory with output alignment folders
cd ../..
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
		#Set output folder name
		outputStats=alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
		#Set header of overall summary csv file
		echo "overall, concordant" > "$outputStats"_allSamples_"$runNum".csv
	elif [[ "$f1" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
		#Set output folder name
		outputStats=alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
		#Set header of overall summary csv file
		echo "mappedLeft, mappedRight, overall, concordant" > "$outputStats"_allSamples_"$runNum".csv
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	#Retrieve summaries for each aligned sample
	for f2 in "$f1"/*/; do
		#Retrieve sample name
		sampleName=$(basename "$f2")
		echo "Merging sample $sampleName alignment summary..."
		#Retrieve sample summary based on alignment method
		bash TranscriptomeAnalysisPipeline_DaphniaUVTolerance/AlignmentAnalysis/alignmentSummary_"$analysisMethod"_sample.sh "$f2" "$analysisMethod" "$runNum"
		#Combine summaries into one csv file
		cat "$outputStats"_combined_"$runNum".csv >> "$outputStats"_allSamples_"$runNum".csv
		rm "$outputStats"_combined_"$runNum".csv
		echo "Sample $sampleName alignment summary has been merged!"
	done
done