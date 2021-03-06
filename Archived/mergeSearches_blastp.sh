#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash mergeSearches_blastp.sh transcriptomeFasta genotype
#Usage Ex: bash mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_proteins
#Usage Ex: bash mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_cds
#Usage Ex: bash mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_transcripts
#Usage Ex: bash mergeSearches_blastp.sh PA42_transcripts PA42_transcripts PA42_proteins
#Usage Ex: bash mergeSearches_blastp.sh PA42_cds PA42_cds PA42_proteins
#Usage Ex: bash mergeSearches_blastp.sh PA42_proteins PA42_proteins PA42_transcripts
#Usage Ex: bash mergeSearches_blastp.sh PA42_proteins PA42_proteins PA42_cds

if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	geno="$2"
	#Set outputs absolute path
	outputFolder="$inputsPath"/"$1"/reciprocalSearched_blastp_"$3"
	#Set blast result paths
	inputDBPath="$outputFolder"/"blastp.outfmt6"
	inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"
elif [[ "$1" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	outputFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
	#Set blast result paths
	inputDBPath="$outputFolder"/"blastp.outfmt6"
	inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	outputFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
	#Set blast result paths
	inputDBPath="$outputFolder"/"blastp.outfmt6"
	inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	outputFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
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
#outFileResults="$outputFolder"/"blastp_merged_summary.txt"
outFileMerged="$outputFolder"/"blastp_merged.txt"
outFileCleaned="$outputFolder"/"blastp_noDuplicates.txt"
cut -f 1,2 "$inputDBPath" > "$outFileCleaned"
awk '{print $2 " " $1}' "$inputRDBPath" > "$outFileCleaned"

#Remove extra tabs
unexpand -a "$outFileCleaned" > "$outFileMerged"

#Output all the duplicates n-1 times
awk 'seen[$0]++' "$outFileMerged" >  "$outFileCleaned"

#Check number of lines
#echo "Recodring number of entries..."
#echo "query,db,queryHits,reciprocalHits,bestHits" > "$outFileResults"
genes1a=$(wc -l "$inputDBPath" | cut -d ' ' -f 1)
genes2a=$(wc -l "$inputRDBPath" | cut -d ' ' -f 1)
genes3b=$(wc -l "$outFileCleaned" | cut -d ' ' -f 1)
echo "$geno","$3","$genes1a","$genes2a","$genes3b"
#echo "Number of entries recorded!"
