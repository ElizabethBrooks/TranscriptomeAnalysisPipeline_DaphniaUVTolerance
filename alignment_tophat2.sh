#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2
#$ -pe smp 1

#Prepare for alignment
cd ..
mkdir aligned_tophat2
module load bio

bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_tophat2/aligned

#Store forward and reverse reads in arrays
FORWARDARRAY=(trimmed/*pForward.fq.gz)
REVERSEARRAY=(trimmed/*pReverse.fq.gz)

#Loop through all forward and reverse paired reads and run tophat2 on each pair
for f1 in trimmed/*pairedForward.fq.gz; do
	tophat2 -G /afs/crc.nd.edu/group/hoth/echo_base/genome/dpulex-genepredict-v11.gff -o aligned_tophat2 aligned_tophat2/aligned ${FORWARDARRAY[*]} ${REVERSEARRAY[*]}
