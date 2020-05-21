#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_genomeBased_transdecoder_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub decoding_genomeBased_transdecoder.sh alignedSequencesFolder
#Usage Ex: qsub decoding_genomeBased_transdecoder.sh aligned_tophat2_run2

#Load necessary modules for ND CRC servers
module load bio/2.0
#Retrieve genome reference and features paths
#genomeRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
multiFASTA=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve TransDecoder software path
softPath=$(grep "transdecoder:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transdecoder://g")
#Retrieve outputs absolute path
inOutputsPath=$(grep "aligning:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
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
../util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     "$genomeFeat" \
     "$multiFASTA" > transcripts.fasta.transdecoder.genome.gff3