#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generate_multiFASTA_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub generate_multiFASTA.sh dataSetName
#Usage Ex: qsub generate_multiFASTA.sh DMelUV
#Note that the data set name input is for output file naming purposes only

#Retrieve genome reference and features paths
pairedReads=$(grep "pairedReads:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pairedReads://g")
#Retrieve outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
outFolder="$outputsPath"
mkdir "$outFolder"
#Combine paired-end read fasta files into a multi fasta
cat "$pairedReads"/*.fq.gz > "$outFolder"/multiFASTA_"$1".fa