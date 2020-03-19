#!/bin/bash
#Script to convert an input gff file to gff3, and expose any issues with the -E flag
#Usage: bash convert_gff_to_gff3.sh genomeVersion
#Usage Ex: bash convert_gff_to_gff3.sh PA42_v3.0

#Load necessary modules
module load bio/cufflinks
#Retrieve input absolute path
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Set output absolute path
outFolder="$genomeFeat"
#Run cufflinks gffread to convert gff file to gff3
gffread -E "$genomeFeat" -o- > "$outFolder"/genomeConverted_"$1".gff3