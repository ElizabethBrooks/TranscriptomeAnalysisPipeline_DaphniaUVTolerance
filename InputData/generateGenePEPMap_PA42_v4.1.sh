#!/bin/bash
#Script to generate protein to gene map
#Usage: bash generateGenePEPMap_PA42_v4.1.sh

ref=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
path=$(dirname $ref)

grep ">" "$path" | sed "s/-pep//g" | sed "s/>//g" > "$path"/col1.txt
tr -d "\n\r" < "$path"/col1.txt | sed "s/dp\_gene/\ndp\_gene/g" > "$path"/col1_cleaned.txt

grep ">" "$path" | sed "s/>//g" > "$path"/col2.txt
tr -d "\n\r" < "$path"/col2.txt | sed "s/dp\_gene/\ndp\_gene/g" > "$path"/col2_cleaned.txt

paste "$path"/col1_cleaned.txt "$path"/col2_cleaned.txt > "$path"/cols_merged.txt

pepMap=$(grep "genePEPMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genePEPMap://g")
tail -n +2 "$path"/cols_merged.txt > "$pepMap"

#Clean up 
rm "$path"/col*.txt