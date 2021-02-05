#!/bin/bash
#Script to convert a genbank fna file to fasta format
#Usage: bash convertFNA_toFASTA.sh fnaFilePath
#Usage ex: bash convertFNA_toFASTA.sh /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fna

#Remove the extension from the input file
fileName=$(echo "$1" | sed 's/\.fna//g')

#First, the fna file is read into std out
#Then, for each sequence the genbank tag and string is removed up to the scaffold tag
#Finally, the string after the scaffold number to the end of the sequence name line is removed
cat "$1" | sed "s/^>.*scaffold/>scaffold/g" | sed "s/,.*//g" > "$fileName".fasta