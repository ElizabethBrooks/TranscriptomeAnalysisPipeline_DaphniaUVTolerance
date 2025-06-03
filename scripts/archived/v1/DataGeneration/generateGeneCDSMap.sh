#!/bin/bash
#Script to create gene to cds map from a cds file
#Usage: generateGeneCDSMap.sh CDSFile genomeTag

inputFile=$1
path=$(dirname $inputFile)

grep ">" $inputFile | sed "s/-*//g" | sed "s/>//g" > $path"/col1.txt"
tr -d "\n\r" < $path"/col1.txt" | sed "s/dp\_gene/\ndp\_gene/g" > $path"/col1_cleaned.txt"

grep ">" $inputFile | sed "s/>//g" > $path"/col2.txt"
tr -d "\n\r" < $path"/col2.txt" | sed "s/dp\_gene/\ndp\_gene/g" > $path"/col2_cleaned.txt"

paste $path"/col1_cleaned.txt" $path"/col2_cleaned.txt" > $path"/cols_merged.txt"

tail -n +2 $path"/cols_merged.txt" > $path"/"$2".gene_cds_map"

#Clean up 
rm $path"/col*.txt"