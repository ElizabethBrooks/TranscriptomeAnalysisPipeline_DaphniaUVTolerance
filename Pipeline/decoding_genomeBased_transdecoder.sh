#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_genomeBased_transdecoder_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub decoding_genomeBased_transdecoder.sh alignedSequencesFolder genomeVersion
#Usage Ex: qsub decoding_genomeBased_transdecoder.sh aligned_tophat2_run2 PA42_v3.0

#Load necessary modules for ND CRC servers
module load bio/2.0
#Retrieve genome reference and features paths
#genomeRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
multiFASTAPath=$(grep "multiFASTA:" ../InputData/inputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
multiFASTA="$multiFASTAPath"/"$1"_multiFASTA.fa
#Retrieve TransDecoder software path
softPath=$(grep "transdecoder:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transdecoder://g")
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outFolder="$outputsPath"/decoded_genomeBased_"$2"
mkdir "$outFolder"
#Move to output folder
cd "$outFolder"
#Clean up genome features file
#sed -e "s/\r//g" "$genomeFeat" > "$outFolder"/genomeFeat_"$2".gff
#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t "$multiFASTA"
#Optionally, identify peptides with homology to known proteins
#TransDecoder.Predict -t transcripts.fasta [ homology options ]