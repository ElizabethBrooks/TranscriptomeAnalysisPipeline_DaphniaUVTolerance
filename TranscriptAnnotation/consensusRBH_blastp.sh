#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash consensusRBH_blastp.sh transcriptomeFasta genotype PA42File outFile
#Usage Ex: bash consensusRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_proteins trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_cds trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh trimmed_run1E05_assemblyTrinity E05 PA42_transcripts trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh PA42_transcripts PA42_transcripts PA42_proteins trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh PA42_cds PA42_cds PA42_proteins trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh PA42_proteins PA42_proteins PA42_transcripts trimmed_run1_PA42_proteins_RBH_summary.txt
#Usage Ex: bash consensusRBH_blastp.sh PA42_proteins PA42_proteins PA42_cds trimmed_run1_PA42_proteins_RBH_summary.txt

#set output summary file path
outFile="$4"

if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input query
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	geno="$2"
	#Set outputs absolute path
	inputQueryFolder="$inputsPath"/"$1"/reciprocalSearched_blastp_"$3"
elif [[ "$1" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	inputQueryFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	inputQueryFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	geno="$2"
	#Set outputs absolute path
	inputQueryFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Determine input DB
if [[ "$3" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	#Set outputs absolute path
	inputDBFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
elif [[ "$3" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	#Set outputs absolute path
	inputDBFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
elif [[ "$3" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	inputsPath=$(dirname "$inputsPath")
	#Set outputs absolute path
	inputDBFolder="$inputsPath"/reciprocalSearched_blastp_"$3"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi

#Merge blast search results
queryFileRBH="$inputQueryFolder"/"blastp_RBH.txt"
dbFileRBH="$inputDBFolder"/"blastp_RBH.txt"

#Set outputs
outFile=$(echo "$outFile" | sed 's/summary/consensusSummary/g')
outFileRBH=$(echo "$outFile" | sed "s/consensusSummary/$geno_consensusRBH/g")

#Pre-clean up
rm $outFileRBH

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" $dbFileRBH; then #RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < $queryFileRBH

#Check number of lines
#echo "query,db,queryRBH,dbRBH,consensusRBH" > "$outFile"
queryHits=$(wc -l "$queryFileRBH" | cut -d ' ' -f 1)
dbHits=$(wc -l "$dbFileRBH" | cut -d ' ' -f 1)
consensusHits=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "$geno","$3","$queryHits","$dbHits","$consensusHits" >> "$outFile"
