#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_genomeBased_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub decoding_transdecoder_genomeBased.sh alignedSequencesFolder
#Usage Ex: qsub decoding_transdecoder_genomeBased.sh aligned_tophat2_run2
#Usage Ex: qsub decoding_transdecoder_genomeBased.sh aligned_hisat2_run1

#Load necessary modules for ND CRC servers
module load bio/2.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine current directory path
currPath=$(dirname $PWD)
#Retrieve genome reference and features paths
#genomeRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
multiFASTA=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve outputs absolute path
inOutputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
outputFolder="$inOutputsPath"/"$1"/decoded_transdecoder_genomeBased
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Clean up genome features file
#sed -e "s/\r//g" "$genomeFeat" > "$outputFolder"/genomeFeat_"$2".gff
#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t "$multiFASTA"
#Optionally, identify peptides with homology to known proteins
TransDecoder.Predict -t "$multiFASTA"
#Generate a genome-based coding region annotation file
"$currPath"/util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     "$genomeFeat" \
     "$multiFASTA" > transcripts.fasta.transdecoder.genome.gff3