#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash consensusRBH.sh proteomeQuery proteomeDB comparisonFile
#Usage ex: bash consensusRBH.sh uvResponsive/Dmel_Svetec_2016 dnaDamageResponse/Dmel_Svetec_2016 PA42_v4.1_proteins

#set output summary file path
dbTag=$(echo "$1" | sed 's/\//_/g')
queryTag=$(echo "$2" | sed 's/\//_/g')

#Set outputs
outPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")
outFile="$outPath"/RBHB/"$dbTag"_"$queryTag"_consensusSummary.txt
outFileRBH="$outPath"/RBHB/"$dbTag"_"$queryTag"_consensusRBH.txt

if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input query
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	dbPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	dbPath="$dbPath"/"$1"
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	dbPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	dbPath="$dbPath"/"$1"
elif [[ "$1" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "proteinSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/proteinSequences://g")
	dbPath=$(dirname "$dbPath")
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "codingSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/codingSequences://g")
	dbPath=$(dirname "$dbPath")
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	dbPath=$(grep "transcriptSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/transcriptSequences://g")
	dbPath=$(dirname "$dbPath")
else
	dbPath=$(grep "databases:" ../InputData/inputPaths.txt | tr -d " " | sed "s/databases://g")
	dbPath="$dbPath"/"$1"
fi

#Determine input DB
if [[ "$2" == *assemblyTrinity* || "$2" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	queryPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	queryPath="$queryPath"/"$2"
elif [[ "$2" == *assembly*Trinity* || "$2" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	queryPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	queryPath="$queryPath"/"$2"
elif [[ "$2" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	queryPath=$(grep "proteinSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/proteinSequences://g")
	queryPath=$(dirname "$queryPath")
elif [[ "$2" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	queryPath=$(grep "codingSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/codingSequences://g")
	queryPath=$(dirname "$queryPath")
elif [[ "$2" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	queryPath=$(grep "transcriptSequences:" ../InputData/inputsPath.txt | tr -d " " | sed "s/transcriptSequences://g")
	queryPath=$(dirname "$queryPath")
else
	#Retrieve database absolute path for querying
	queryPath=$(grep "databases:" ../InputData/inputPaths.txt | tr -d " " | sed "s/databases://g")
	queryPath="$queryPath"/"$2"
fi

#Merge blast search results
consensusTag=$(echo "$3" | sed 's/\//_/g')
queryFileRBH="$queryPath"/reciprocalSearched_blastp_"$consensusTag"/"blastp_RBH.txt"
dbFileRBH="$dbPath"/reciprocalSearched_blastp_"$consensusTag"/"blastp_RBH.txt"

#Report inputs
echo "Query: $queryFileRBH"
echo "DB: $dbFileRBH"

#Remove headers
tail -n +2 $queryFileRBH > tmp1.txt
tail -n +2 $dbFileRBH > tmp2.txt

#Pre-clean up
echo "queryHit,dbHit" > $outFileRBH

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q ",$f2" tmp2.txt; then #consensus RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < tmp1.txt

#Check number of lines
echo "query,db,consesus,queryRBH,dbRBH,consensusRBH" > "$outFile"
queryHits=$(($(wc -l "$queryFileRBH" | cut -d ' ' -f 1)-1))
dbHits=$(($(wc -l "$dbFileRBH" | cut -d ' ' -f 1)-1))
consensusHits=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "$1","$2","$3","$queryHits","$dbHits","$consensusHits" >> "$outFile"

#Clean up
rm tmp*.txt