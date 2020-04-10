#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_overall.sh alignmentFolders
#Usage Ex: bash alignmentSummary_overall.sh aligned_hisat2_run2 aligned_tophat2_run1 aligned_tophat2_run3
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve input alignment summary absolute path
inputsPath=$(grep "aligning:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory
outDir="$outputsPath"/AlignmentsAnalyzed
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
		outputStats="$outDir"/alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
		#Set header of overall summary csv file
		echo "sample,overall,concordant" > "$outputStats"_"$runNum".csv
	elif [[ "$f1" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
		#Set output folder name
		outputStats="$outDir"/alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_//g")
		#Set header of overall summary csv file
		echo "sample,leftMapped,rightMapped,overall,concordant" > "$outputStats"_"$runNum".csv
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	echo "Merging $f1 alignment summaries..."
	#Retrieve summaries for each aligned sample
	for f2 in "$inputsPath"/"$f1"/*; do
		#Retrieve sample name
		sampleName=$(basename "$f2")
		#Retrieve sample summary based on alignment method
		bash alignmentSummary_"$analysisMethod"_sample.sh "$f2" "$analysisMethod" "$runNum"
		#Combine summaries into one csv file
		cat "$outputStats"_combined_"$runNum".csv >> "$outputStats"_"$runNum".csv
		rm "$outputStats"_combined_"$runNum".csv
	done
	echo "Alignment summaries for $f1 have been merged!"
done