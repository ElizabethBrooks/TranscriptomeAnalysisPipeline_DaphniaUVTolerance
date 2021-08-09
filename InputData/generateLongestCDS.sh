#!/bin/bash
#Script to generate the longest CDS from the PA42 v4.1 gff and reference fasta
#Usage: bash generateLongestCDS.sh

#Load necessary module
module load bio

#Retrieve input files
inputFeat=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
inputRef=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Set output file names
outCDS=$(dirname "$inputRef")
outCDS="$outCDS"/PA42_v4.1_longestCDS.fa

#Usage: gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] [-o <outfile.gff>] [-t <tname>] [-r #[[<strand>]<chr>:]<start>..<end> [-R]] [-CTVN‐ JMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>] [-i <maxintron>]
gffread "$inputFeat" -g "$inputRef" -x "$outCDS" -W -F
