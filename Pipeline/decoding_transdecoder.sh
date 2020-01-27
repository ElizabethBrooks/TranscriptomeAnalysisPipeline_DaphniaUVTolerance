#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to predict coding regions from a transcript fasta file

#Move up a directory
cd ..

#Construct the transcript fasta file using the genome and the transcripts.gtf file
util/gtf_genome_to_cdna_fasta.pl transcripts.gtf test.genome.fasta > transcripts.fasta

#Convert the transcript structure GTF file to an alignment-GFF3 formatted file
util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3

#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t transcripts.fasta
#Optionally, identify peptides with homology to known proteins
TransDecoder.Predict -t transcripts.fasta [ homology options ]

#Generate a genome-based coding region annotation file
util/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
