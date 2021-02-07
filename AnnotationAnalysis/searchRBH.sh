#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchRBH_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: qsub searchRBH.sh transcriptomeFasta genotype PA42File outFile
#Usage Ex: qsub searchRBH.sh trimmed_run1E05_assemblyTrinity E05 PA42_v3.0_proteins rbhb_summary.txt
#Usage Ex: qsub searchRBH.sh trimmed_run1E05_assemblyTrinity E05 PA42_v3.0_cds rbhb_summary.txt
#Usage Ex: qsub searchRBH.sh trimmed_run1E05_assemblyTrinity E05 PA42_v3.0_transcripts rbhb_summary.txt
#Usage Ex: qsub searchRBH.sh PA42_v4.1_transcripts PA42_v4.1_transcripts PA42_v4.1_proteins rbhb_consensusSummary.txt rbhb_uniqueRBH.txt
#Usage Ex: qsub searchRBH.sh PA42_v4.1_cds PA42_cds PA42_v4.1_proteins rbhb_consensusSummary.txt rbhb_uniqueRBH.txt
#Usage Ex: qsub searchRBH.sh PA42_v4.1_proteins PA42_v4.1_proteins PA42_v4.1_transcripts rbhb_consensusSummary.txt rbhb_uniqueRBH.txt
#Usage Ex: qsub searchRBH.sh PA42_v4.1_proteins PA42_v4.1_proteins PA42_v4.1_cds rbhb_consensusSummary.txt rbhb_uniqueRBH.txt

if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input database for blastp
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	geno="$2"
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	geno="$2"
elif [[ "$1" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$1"
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$1"
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$1"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Set outputs absolute path
outputFolder="$inputsPath"/"$1"/reciprocalSearched_blastp_"$3"
#Set blast result paths
inputDBPath="$outputFolder"/"blastp.outfmt6"
inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"

#Move to output folder
cd "$outputFolder"

#set output paths
outFile="$4"
outFileRBH="$outputFolder"/"blastp_RBH.txt"

#Report inputs
echo "Query: $inputDBPath"
echo "DB: $inputRDBPath"

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
