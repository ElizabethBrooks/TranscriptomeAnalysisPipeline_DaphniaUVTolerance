#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_genomeGuided_trinity_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub assembly_genomeGuided_trinity.sh
#Usage Ex: qsub assembly_genomeGuided_trinity.sh
#Genome-guided Trinity De novo Transcriptome Assembly

#TO DO: double check required modules
#module load bio
module load bio/trinity/2.8.4
#module load bio/salmon/0.11.2
#module load bio/jellyfish/2.2.10
#module load bio/python/2.7.14

#Users must provide read alignments to Trinity as a coordinate-sorted bam file
#Be sure it's coordinate sorted by running 'samtools sort' on it
#Use a maximum intron length that makes most sense given your targeted organism
Trinity --genome_guided_bam rnaseq.coordSorted.bam \
         --genome_guided_max_intron 10000 \
         --max_memory 10G --CPU 10 