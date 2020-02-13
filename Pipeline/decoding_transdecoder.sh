#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: bash decoding_transdecoder.sh genomeVersion
#Usage Ex: bash decoding_transdecoder.sh PA42_v3.0

#Load necessary modules for CRC servers
module load bio/transdecoder/
#Retrieve genome reference path
genomeRef=$(grep "genomeReference:" ../InputData/outputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFeat=$(grep "genomeFeatures:" ../InputData/outputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Construct the transcript fasta file using the genome and the transcripts.gtf file
../util/gtf_genome_to_cdna_fasta.pl "$genomeFeat" "$genomeRef" > transcripts_"$2".fasta
#Convert the transcript structure GTF file to an alignment-GFF3 formatted file
../util/gtf_to_alignment_gff3.pl "$genomeFeat" > transcripts_"$2".gff3
#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t transcripts_"$2".fasta
#Optionally, identify peptides with homology to known proteins
#TransDecoder.Predict -t transcripts.fasta [ homology options ]
#Generate a genome-based coding region annotation file
../util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder_"$2".gff3 \
     transcripts_"$2".gff3 \
     transcripts_"$2".fasta > transcripts.fasta.transdecoder.genome_"$2".gff3
