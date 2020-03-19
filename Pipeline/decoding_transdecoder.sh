#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#Script to predict coding regions from a transcript fasta file using Transdecoder
#Usage: qsub decoding_transdecoder.sh genomeVersion
#Usage Ex: qsub decoding_transdecoder.sh PA42_v3.0
#Note that the genome version input is for output file naming purposes only

#Load necessary modules for ND CRC servers
module load bio/transdecoder
#module load bio/cufflinks
#Retrieve genome reference and features paths
multiFASTAPath=$(grep "multiFASTA:" ../InputData/inputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
multiFASTA="$multiFASTAPath"/multiFASTA_"$1".fa
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outFolder="$outputsPath"/decoded_"$1"
mkdir "$outFolder"
#Move to output folder
cd "$outFolder"
#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t "$multiFASTA"
#Optionally, identify peptides with homology to known proteins
#TransDecoder.Predict -t transcripts.fasta [ homology options ]