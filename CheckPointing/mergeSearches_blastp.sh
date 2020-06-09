#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash mergeSearches_blastp.sh transcriptomeFasta genotype
#Usage Ex: bash mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searchedTranscripts_blastp
	#Set blast result paths
	inputDBPath="$outputFolder"/"blastp.outfmt6"
	inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Move to output folder
cd "$outputFolder"

#Merge blast search results
outFileResults="$outputFolder"/"blastp_merged_summary.txt"
outFileMerged="$outputFolder"/"blastp_merged.txt"
outFileCleaned="$outputFolder"/"blastp_noDuplicates.txt"
cut -f 1,2 "$inputDBPath" > "$outFileCleaned"
awk '{print $2 " " $1}' "$inputRDBPath" > "$outFileCleaned"

#Remove extra tabs
unexpand -a "$outFileCleaned" > "$outFileMerged"

#Output all the duplicates n-1 times
awk 'seen[$0]++' "$outFileMerged" >  "$outFileCleaned"

#Check number of lines
echo "Recodring number of entries..."
echo "Query, Total, Unique, Duplicates" > "$outFileResults"
genes1a=$(wc -l "$inputDBPath" | cut -d ' ' -f 1)
genes1b=$(sort "$inputDBPath" | uniq -u | wc -l)
genes1c=$(($genes1a-$genes1b))
echo "$2"", $genes1a, $genes1b, $genes1c" >> "$outFileResults"
genes2a=$(wc -l "$inputRDBPath" | cut -d ' ' -f 1)
genes2b=$(sort "$inputRDBPath" | uniq -u | wc -l)
genes2c=$(($genes2a-$genes2b))
echo "Pulex, $genes2a, $genes2b, $genes2c" >> "$outFileResults"
genes3a=$(wc -l "$outFileMerged" | cut -d ' ' -f 1)
genes3b=$(wc -l "$outFileCleaned" | cut -d ' ' -f 1)
genes3c=$(($genes3a-$genes3b))
echo "Merged, $genes3a, $genes3c, $genes3b" >> "$outFileResults"
echo "Number of entries recorded!"
