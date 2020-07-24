#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N mergeSearches_blastp_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: qsub mergeSearches_blastp.sh transcriptomeFasta genotype PA42File outFile
#Usage Ex: qsub mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_cds trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_transcripts trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh PA42_transcripts PA42_transcripts PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh PA42_cds PA42_cds PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh PA42_proteins PA42_proteins PA42_transcripts trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub mergeSearches_blastp.sh PA42_proteins PA42_proteins PA42_cds trimmed_run1_PA42_proteins_blastp_summary.txt

#set output summary file path
outFile="$4"

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
outFileRBH="$outputFolder"/"blastp_RBH.txt"
outFileQuery="$outputFolder"/"tmp1_blastp.txt"
outFileDB="$outputFolder"/"tmp2_blastp.txt"
awk '{print $1 "," $2}' "$inputDBPath" > "$outFileQuery"
awk '{print $2 "," $1}' "$inputRDBPath" > "$outFileDB"

#Pre-clean up
echo "query,db" > $outFileRBH

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" $outFileDB; then #RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < $outFileQuery

#Check number of lines
#echo "Recodring number of entries..."
#echo "query,db,queryHits,dbHits,bestHits" > "$outFileResults"
genes1a=$(wc -l "$outFileQuery" | cut -d ' ' -f 1)
genes2a=$(wc -l "$outFileDB" | cut -d ' ' -f 1)
genes3b=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "$geno","$3","$genes1a","$genes2a","$genes3b" >> "$outFile"

#Clean up
rm "$outputFolder"/tmp*
