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
#Prepare tags file for comparison
grep "Pool_1" mergeGuideFile_"$1"_"$2".txt > tmp.txt | sed -ie 's|^|Pool_1_|' tmp.txt
grep "Pool_2" mergeGuideFile_"$1"_"$2".txt >> tmp.txt | sed -ie 's|^|Pool_2_|' tmp.txt
grep "Pool_3" mergeGuideFile_"$1"_"$2".txt >> tmp.txt | sed -ie 's|^|Pool_3_|' tmp.txt
#Loop through all counted paired reads and append each sample tag
# with the corresponding file path
cat ../InputData/mergeCounts_tags_"$2".txt > "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
for f1 in "$inputsPath"/"$1"/*/; do
	currSample=$(basename "$f1" | sed "s/140327_I481_FCC3P1PACXX_L..//g")
	currTag=$(grep "$currSample" tmp.txt | sed "s/^Pool_.//g")
	sed -i 's/'"$currTag"'/'"$currTag"' '"$f1"'/g' "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
done
#Merge gene counts based on generated guide file
python merge_tables.py "$outputsPath"/mergeGuideFile_"$1"_"$2".txt
#Rename the output merged counts file
mv merged_counts.txt "$outputsPath"/mergedCounts_"$1"_"$2".txt
