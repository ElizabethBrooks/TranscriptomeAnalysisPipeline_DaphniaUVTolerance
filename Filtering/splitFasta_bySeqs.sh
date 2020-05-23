#!/bin/bash
#Script to split an input fasta file input appropriate number of sequence
# chunks based on the input max chunk size (in MB)
#Usage: bash splitFasta_bySeqs.sh fastaFile maxMB
#Usage Ex: bash splitFasta_bySeqs.sh ../../Daphnia_pulex_PA42_v3.0.fasta 50

#Retrieve input fasta absolute path
inputsPath=$(dirname "$1")
#Determine the number of pieces to split the file into
chunkSize=$2 #Max number of MB for each chunk
#Determine the size of the input fasta file
fastaSize=$(ls -l --block-size=1MB intronDataDriver_fromGFF3.sh | cut -d " " -f5)
chunkNum=$((($fastaSize+($chunkSize-1))/$chunkSize)) #Round up
#Determine the number of sequences to include in each chunk
seqsNum=$(grep ">" "$1" | wc -l)
seqsNum=$((($seqsNum+($chunkNum-1))/$chunkNum)) #Round up
#Retrieve the input fatsa file name
inFileNoEx=$(basename "$1")
inFileNoEx=$(echo $inFileNoEx | sed 's/\.fasta//')
#Set output file name
outFile="$inputsPath"/"$inFileNoEx"
#Loop through the input fasta file and split by sequence chunks
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%$seqsNum==0){file=sprintf("myseq%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < "$1"
