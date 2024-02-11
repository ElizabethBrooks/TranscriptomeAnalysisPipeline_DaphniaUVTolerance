#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchRBH_jobOutput
#$ -q debug

# script to filter reciprocal blast results for best hits
# Usage: qsub searchRBH.sh searchTag
# Usage ex: qsub searchRBH.sh Dpulex_Dmelanogaster

if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi

# retrieve inputs absolute path
inputsPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")

# set inputs and outputs absolute folder name
searchTag="$1"
inputsPath=$inputsPath"/reciprocalSearched_blastp_"$searchTag
outputFolder=$inputsPath

#Set blast result paths
inputDBPath="$outputFolder"/"blastp.outfmt6"
inputRDBPath="$outputFolder"/"blastp_reciprocal.outfmt6"

#Move to output folder
cd "$outputFolder"

#set output paths
outFile=$outputFolder"/blastp_consensusSummary.txt"
outFileRBH=$outputFolder"/blastp_RBH.txt"

#Report inputs
echo "Query: $inputDBPath"
echo "DB: $inputRDBPath"

#Pre-clean up to find best hits
bestHitsDB="$outputFolder"/"blastp_bestHits.outfmt6"
bestHitsRDB="$outputFolder"/"blastp_bestHits_reciprocal.outfmt6"
awk '!seen[$1]++' $inputDBPath > $bestHitsDB
awk '!seen[$1]++' $inputRDBPath > $bestHitsRDB

#Pre-clean up
echo "queryHit,dbHit,db" > $outFileRBH

#Loop over first set of annotations
while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12
do
	#Determine annotation sets
	if grep -q "$f2"$'\t'"$f1"$'\t' $bestHitsRDB; then #RBH
		echo "$f1,$f2" >> $outFileRBH
	fi
done < $bestHitsDB

#Check number of lines
#echo "Recodring number of entries..."
#echo "query,db,queryHits,dbHits,bestHits,similarity" > "$outFileResults"
queryHits=$(wc -l "$bestHitsDB" | cut -d ' ' -f 1)
dbHits=$(wc -l "$bestHitsRDB" | cut -d ' ' -f 1)
bestHits=$(($(wc -l "$outFileRBH" | cut -d ' ' -f 1)-1))
echo "$2","$3","$queryHits","$dbHits","$bestHits" >> "$outFile"
