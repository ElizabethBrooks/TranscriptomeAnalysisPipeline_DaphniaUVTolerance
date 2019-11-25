#!/bin/bash
#Usage: bash geneCounts_generateMergedTable.sh countedGenesFolder sampleSet
#Usage Ex: bash geneCounts_generateMergedTable.sh counted_htseqTophat2_run1 fullset
#Script to generate guide file and merge gene counts
# using the merge_tables.py script
#Retrieve path for input counted reads
inputsPath=$(grep "counting:" InputData/outputPaths.txt | tr -d " " | sed "s/counting://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Loop through all counted paired reads and append each sample tag
# with the corresponding file path
for f1 in "$1"/*/; do
	sed -i -e 's|^|'"$inputsPath"'/'"$f1"'counts.txt|' "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
done
python merge_tables.py "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
#Rename the output merged counts file
mv merged_counts.txt "$outputsPath"/mergedCounts_"$1"_"$2".txt
