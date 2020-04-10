#!/bin/bash
#Usage: bash geneCounts_generateMergedTable.sh countedGenesFolder sampleSet
#Usage Ex: bash geneCounts_generateMergedTable.sh counted_htseqTophat2_run1 fullset
#Script to generate guide file and merge gene counts
# using the merge_tables.py script
#Load necessary modules for ND CRC servers
module load bio/python/2.7.14
#Prepare for analysis
dirFlag=0
runNum=1
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
inputsPath=$(grep "counting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/counting://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Make a new directory for the analysis
outputsPath="$outputsPath"_"$1"_"$2"/GeneCountsAnalyzed
mkdir $outputsPath
echo "Creating folder for gene count analysis for a $2 of $1..."
#Remove extra number tags from file names
for f0 in "$inputsPath"/"$1"/*/; do
	newName=$(echo "$f0" | sed 's/UV1/UV/g')
	newName=$(echo "$newName" | sed 's/VIS1/VIS/g')
	newName=$(echo "$newName" | sed 's/UV2/UV/g')
	newName=$(echo "$newName" | sed 's/VIS2/VIS/g')
	newName=$(echo "$newName" | sed 's/UV3/UV/g')
	newName=$(echo "$newName" | sed 's/VIS3/VIS/g')
	if [[ "$newName" != "$f0" ]]; then
		mv "$f0" "$newName"
	fi
done
#Prepare tags file for comparison
grep "Pool1" ../InputData/mergeCounts_guideFile_tags_"$2".txt | sed 's/^/Pool_1_/' > "$outputsPath"/tmp1.txt
grep "Pool2" ../InputData/mergeCounts_guideFile_tags_"$2".txt | sed 's/^/Pool_2_/' > "$outputsPath"/tmp2.txt
grep "Pool3" ../InputData/mergeCounts_guideFile_tags_"$2".txt | sed 's/^/Pool_3_/' > "$outputsPath"/tmp3.txt
cat "$outputsPath"/tmp*.txt >> "$outputsPath"/tmp.txt
#Loop through all counted paired reads and append each sample tag
# with the corresponding file path
cat ../InputData/mergeCounts_guideFile_tags_"$2".txt > "$outputsPath"/mergeCounts_guideFile_"$1"_"$2".txt
for f1 in "$inputsPath"/"$1"/*/; do
	currSample=$(basename "$f1" | sed "s/140327_I481_FCC3P1PACXX_L..//g")
	#Determine if subset of files are to be used
	if grep -iFq "$currSample" "$outputsPath"/tmp.txt; then
		currTag=$(grep "$currSample" "$outputsPath"/tmp.txt | sed "s/Pool_._//g")
		sed -i 's,'"$currTag"','"$f1"'counts.txt '"$currTag"',' "$outputsPath"/mergeCounts_guideFile_"$1"_"$2".txt
	fi
done
#Clean up
rm "$outputsPath"/tmp*.txt
#Move to location of merge_tagles.py script
cd ../util
#Merge gene counts based on generated guide file
python merge_tables.py "$outputsPath"/mergeCounts_guideFile_"$1"_"$2".txt
#Rename the output merged counts file
mv merged_counts.txt "$outputsPath"/geneCounts_merged_"$1"_"$2".txt
#Print a script completion confirmation message
echo "Merged table has been renamed 'geneCounts_merged_"$1"_"$2".txt' and moved!"