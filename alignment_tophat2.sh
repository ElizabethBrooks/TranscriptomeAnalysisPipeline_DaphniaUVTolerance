#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2_jobOutput
#$ -pe smp 1

#Prepare for alignment
cd ..
mkdir aligned_tophat2
mkdir aligned_tophat2/out
module load bio
#Build reference genome
bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_tophat2/Daphnia_pulex.allmasked
#Loop through all forward and reverse paired reads and run tophat2 on each pair
for f1 in trimmed/*pForward.fq.gz; do
	tophat2 -G /afs/crc.nd.edu/group/hoth/echo_base/genome/dpulex-genepredict-v11.gff -o aligned_tophat2/out/"${f1:8:${#f1}-23}" aligned_tophat2/Daphnia_pulex.allmasked $f1 "${f1:0:${#f1}-14}"pReverse.fq.gz
done
