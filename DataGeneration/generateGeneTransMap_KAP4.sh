#!/bin/bash

# script to create a gene to transcript map for the KAP4 genome using the transcript sequences file
# that have gene ID annotations
# usage ex: bash generateGeneTransMap_KAP4.sh

# retrieve input path
transcriptPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")

# retrieve output map path
mapPath=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")

# set outputs path
path=$(dirname $mapPath)

# retrieve transcript IDs
cat $transcriptPath | grep ">" | cut -d " " -f 1 | sed "s/\>//g" > $path"/col1_cleaned.txt"

# retrieve gene IDs
cat $transcriptPath | grep ">" | cut -d "(" -f 2 | cut -d ")" -f 1 > $path"/col2_cleaned.txt"

# create txt file with both sets of IDs as columns
paste $path"/col1_cleaned.txt" $path"/col2_cleaned.txt" > $mapPath

# clean up 
rm $path"/col"*"_cleaned.txt"