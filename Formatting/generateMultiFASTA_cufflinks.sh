#!/bin/bash
#Script to generate a multi FASTA file for all transcripts in a GFF file
#Usage: qsub decoding_transdecoder.sh genomeVersion
#Usage Ex: qsub decoding_transdecoder.sh PA42_v3.0
#Note that the genome version input is for output file naming purposes only

#Load necessary modules for ND CRC servers
module load bio/samtools
module load bio/cufflinks
#Retrieve reads and genome features paths
pairedReads=$(grep "pairedReads:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pairedReads://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
outFolder="$outputsPath"
mkdir "$outFolder"
#Generate a fasta index file using samtools
samtools faidx "$outFolder"/genomeIndex_"$1".fa
#Generate a FASTA file with the DNA sequences for all transcripts in the GFF file
gffread -w transcripts.fa -g "$outFolder"/genomeIndex_"$1".fa "$genomeFeat"