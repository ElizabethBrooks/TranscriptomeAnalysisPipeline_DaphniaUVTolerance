#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_test
#$ -pe smp 8

cd ..
mkdir aligned_hisat2_test

#Prepare for mapping
module load bio/hisat2/2.1.0
hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_hisat2_test/Daphnia_pulex.allmasked

#Loop through all forward and reverse paired reads and run hisat2 on each pair
for f1 in trimmed/*pForward.fq.gz; do
	hisat2 -q -p 8 -x aligned_hisat2_test/Daphnia_pulex.allmasked -1 $f1 -2 "${f1:0:${#f1}-14}"pReverse.fq.gz -S aligned_hisat2_test/out/"${f1:0:${#f1}-14}" --summary-file alignedSummary_test.txt
done