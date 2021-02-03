#!/bin/bash
#Script to perform reformatting of normalized count sets to GCT format
#Usage: bash intronDataDriver_fromGFF3.sh
#Usage Ex: bash intronDataDriver_fromGFF3.sh

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Set paths for r script
export R_LIBS=/afs/crc.nd.edu/user/e/ebrooks5/R/x86_64-pc-linux-gnu-library/3.5:$R_LIBS
#Retrieve analysis inputs path
inputsPath=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Generate intron data using GenomicFeatures in R
Rscript generateIntronData_fromGFF3.r "$inputsPath"