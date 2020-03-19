#!/bin/bash
#Script to generate a multi FASTA file from aligned transcript sequences
#Usage: bash generateMultiFASTA_bash.sh alignedSequencesFolder
#Usage Ex: bash generateMultiFASTA_bash.sh aligned_tophat2_run2

#Retrieve genome reference and features paths
pairedReads=$(grep "alignment:" ../InputData/inputPaths.txt | tr -d " " | sed "s/alignment://g")
#Retrieve outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
outFolder="$outputsPath"
mkdir "$outFolder"
#Combine paired-end read fasta files into a multi fasta
cat "$pairedReads"/*.fq.gz > "$outFolder"/"$1"_multiFASTA.fa