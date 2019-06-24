#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 8

#Prepare for mapping
cd ..
mkdir aligned_hisat2
mkdir aligned_hisat2/out
module load bio/hisat2/2.1.0
#Build reference genome
#... add check ...
hisat1-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_hisat2/Daphnia_pulex.allmasked
#Loop through all forward and reverse paired reads and run hisat2 on each pair
# using 8 threads
for f1 in trimmed/*pForward.fq.gz; do
	hisat2 -q -p 8 -x aligned_hisat2/Daphnia_pulex.allmasked -1 $f1 -2 "${f1:0:${#f1}-14}"pReverse.fq.gz -S aligned_hisat2/out/"${f1:8:${#f1}-23}" --summary-file aligned_hisat2/alignedSummary.txt
done
