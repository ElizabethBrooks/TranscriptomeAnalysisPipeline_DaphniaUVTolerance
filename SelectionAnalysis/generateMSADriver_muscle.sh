#!/bin/bash

#Retrieve inputs path
refPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
#Set outputs path
outPath=$(grep "MSA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/MSA://g")
outPath="$outPath"/daphniaMSA_PA42_v4.1
#Check if the folder already exists
mkdir "$outPath"
if [ $? -ne 0 ]; then
	echo "The $outPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Loop over all genes in the reference
grep ">" "$refPath" | sed "s/-CDS//g" | sed "s/>//g" > "$outPath"/col1.txt

#Split the reference gene set into segments
split -l 1000 "$outPath"/col1.txt "$outPath"/colRef

#Generate MSAs for each segment
for i in "$outPath"/colRef*; do
	qsub generateMSA_muscle.sh "$i" "$@"
done

#Clean up
rm "$outPath"/col1.txt