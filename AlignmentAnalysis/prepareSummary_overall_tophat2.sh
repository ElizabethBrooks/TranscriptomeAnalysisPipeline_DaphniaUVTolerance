#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash prepareSummary_overall_tophat2.sh alignmentOutputFolder alignmentMethod
#Usage Ex: bash prepareSummary_overall_tophat2.sh aligned_topaht2_run2 topaht2
#Determine if the folder name was input in the correct format
if [[ $1 == *\/* ]] || [[ $1 == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Prepare input and output file names
inputStats="$1"
outputStats=alignmentStats_"$2"
#Retrieve summaries for each aligned sample
for f1 in "$1"/*/; do
	bash prepareSummary_sample_tophat2.sh "$f1" $2
	#Combine summaries into one csv file
	cat "$outputStats"_combined.csv >> "$outputStats"_allSamples.csv
done