#!/bin/bash
#Usage: bash geneCounts_generateMergedTable.sh countedGenesFolder sampleSet
#Usage Ex: bash geneCounts_generateMergedTable.sh counted_htseqTophat2_run1 fullset
#Script to generate guide file and merge gene counts
# using the merge_tables.py script
#Load necessary modules for ND CRC servers
module load bio/python/2.7.14
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Retrieve path for input counted reads
inputsPath=$(grep "counting:" InputData/outputPaths.txt | tr -d " " | sed "s/counting://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Loop through all counted paired reads and append each sample tag
# with the corresponding file path
for f1 in "$inputsPath"/"$1"/*/; do
	sed -i -e 's|^|'"$f1"'counts.txt|' "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
done
python merge_tables.py "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
#Rename the output merged counts file
mv merged_counts.txt "$outputsPath"/mergedCounts_"$1"_"$2".txt
