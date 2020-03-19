#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_genomeBased_transdecoder_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub decoding_genomeBased_transdecoder.sh genomeVersion
#Usage Ex: qsub decoding_genomeBased_transdecoder.sh PA42_v3.0
#Note that the genome version input is for output file naming purposes only

#Load necessary modules for ND CRC servers
module load bio/transdecoder
module load bio/cufflinks
#Retrieve genome reference and features paths
genomeRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve TransDecoder software path
softPath=$(grep "transdecoder:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transdecoder://g")
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outFolder="$outputsPath"/decoded_genomeBased_"$1"
mkdir "$outFolder"
#Clean up genome features file
#sed -e "s/\r//g" "$genomeFeat" > "$outFolder"/genomeFeat_"$2".gff
#Convert gff to gff3, and expose any issues with the -E flag
gffread -E "$genomeFeat" -o- "$outFolder"/genomeFeat_"$2".gff3
#Move to transdecoder software folder
#cd "$softPath"
#Construct the transcript fasta file using the genome and the transcripts.gtf file (cufflinks of stringtie)
#perl util/gtf_genome_to_cdna_fasta.pl "$outFolder"/genomeFeat_"$2".gtf "$genomeRef" > "$outFolder"/transcripts_"$1".fasta
#Convert the transcript structure GTF file to an alignment-GFF3 formatted file
#perl util/gtf_to_alignment_gff3.pl "$outFolder"/genomeFeat_"$2".gtf > "$outFolder"/transcripts_"$1".gff3
#Generate your best candidate open rading frame (ORF) predictions
#TransDecoder.LongOrfs -t "$outFolder"/transcripts_"$1".fasta
#Optionally, identify peptides with homology to known proteins
#TransDecoder.Predict -t transcripts.fasta [ homology options ]
#Generate a genome-based coding region annotation file
#perl util/cdna_alignment_orf_to_genome_orf.pl \
#     "$outFolder"/transcripts_"$1".fasta.transdecoder.gff3 \
#     "$outFolder"/transcripts_"$1".gff3 \
#     "$outFolder"/transcripts_"$1".fasta > "$outFolder"/transcripts.fasta.transdecoder.genome_"$1".gff3