#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchRBH_blastp_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: qsub searchRBH_blastp.sh transcriptomeFasta genotype PA42File outFile
#Usage Ex: qsub searchRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_cds trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_transcripts trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh PA42_transcripts PA42_transcripts PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh PA42_cds PA42_cds PA42_proteins trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh PA42_proteins PA42_proteins PA42_transcripts trimmed_run1_PA42_proteins_blastp_summary.txt
#Usage Ex: qsub searchRBH_blastp.sh PA42_proteins PA42_proteins PA42_cds trimmed_run1_PA42_proteins_blastp_summary.txt

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

#set output paths
outFile="$4"
outFileRBH="$outputFolder"/"blastp_RBH.txt"


#Pre-clean up
echo "queryHit,dbHit,db" > $outFileRBH

#Loop over first set of annotations
while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12
do
	#Determine annotation sets
	if grep -q "$f2"$'\t'"$f1"$'\t' $inputRDBPath; then #RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < $inputDBPath

#Check number of lines
#echo "Recodring number of entries..."
#echo "query,db,queryHits,dbHits,bestHits,similarity" > "$outFileResults"
queryHits=$(wc -l "$inputDBPath" | cut -d ' ' -f 1)
dbHits=$(wc -l "$inputRDBPath" | cut -d ' ' -f 1)
bestHits=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "$geno","$3","$queryHits","$dbHits","$bestHits" >> "$outFile"
