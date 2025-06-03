#!/bin/bash
#Script to convert an input gff file to gff3, and expose any issues with the -E flag
#Usage: bash convertGFF_toGTF.sh genomeVersion
#Usage Ex: bash convertGFF_toGTF.sh PA42_v3.0

#Load necessary modules
module load bio/cufflinks
#Retrieve input absolute path
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Re set output extension
outFile=$(echo $genomeFeat | sed 's/\.gff/\.gtf/g')
#Run cufflinks gffread to convert gff file to gff3
gffread -E "$genomeFeat" -o- > "$outFile"