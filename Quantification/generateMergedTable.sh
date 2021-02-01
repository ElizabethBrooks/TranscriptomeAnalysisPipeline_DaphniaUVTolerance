#!/bin/bash
#Usage: bash generateMergedTable.sh countedGenesFolder sampleSet
#Usage Ex: bash generateMergedTable.sh counted_htseq_run1 sortedName_samtoolsHisat2_run1 trimmed_run1E05_assemblyTrinity
#Usage Ex: bash generateMergedTable.sh counted_htseq_run1 sortedName_samtoolsHisat2_run2 genome
#Script to generate guide file and merge gene counts
# using the merge_tables.py script
#Load necessary modules for ND CRC servers
module load bio
#module load bio/python/2.7.14
#Prepare for analysis
dirFlag=0
runNum=1
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$3" == *assemblyTrinity* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsPath="$inputsPath"/"$3"/"$2"/"$1"
elif [[ "$1" == *assemblyGenome* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsPath="$inputsPath"/"$3"/"$2"/"$1"
elif [[ "$3" == "genome" ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsPath="$inputsPath"/"$2"/"$1"
else
	echo "ERROR: The counted files were not found... exiting"
	exit 1
fi
#Make a new directory for the analysis
outputsPath="$inputsPath"/GeneCountsAnalyzed
mkdir $outputsPath
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi
#Remove any extra number tags from file names
for f0 in "$inputsPath"/*/; do
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
grep "Pool1" ../InputData/mergeCounts_guideFile_tags.txt | sed 's/^/Pool_1_/' > "$outputsPath"/tmp1.txt
grep "Pool2" ../InputData/mergeCounts_guideFile_tags.txt | sed 's/^/Pool_2_/' > "$outputsPath"/tmp2.txt
grep "Pool3" ../InputData/mergeCounts_guideFile_tags.txt | sed 's/^/Pool_3_/' > "$outputsPath"/tmp3.txt
cat "$outputsPath"/tmp*.txt >> "$outputsPath"/tmp.txt
#Loop through all counted paired reads and append each sample tag
# with the corresponding file path
guideFile="$outputsPath"/tmp_mergeCounts_guideFile_"$3"_"$2"_"$1".txt
cat ../InputData/mergeCounts_guideFile_tags.txt > "$guideFile"
for f1 in "$inputsPath"/*/; do
	currSample=$(basename "$f1" | sed "s/140327_I481_FCC3P1PACXX_L..//g")
	#Determine if subset of files are to be used
	if grep -iFq "$currSample" "$outputsPath"/tmp.txt; then
		currTag=$(grep "$currSample" "$outputsPath"/tmp.txt | sed "s/Pool_._//g")
		sed -i 's,'"$currTag"','"$f1"'counts.txt '"$currTag"',' "$guideFile"
	fi
done
#Move to location of merge_tagles.py script
cd ../util
#Merge gene counts based on generated guide file
python merge_tables.py "$guideFile"
#Rename the output merged counts file
mergedCounts="$outputsPath"/geneCounts_merged_"$3"_"$2"_"$1".txt
mv merged_counts.txt "$mergedCounts"
#Print a script completion confirmation message
echo "Merged table has been renamed $mergedCounts and moved!"
#Clean up
rm "$outputsPath"/tmp*.txt