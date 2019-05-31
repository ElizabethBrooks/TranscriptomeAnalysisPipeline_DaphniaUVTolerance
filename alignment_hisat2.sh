#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2
#$ -pe smp 1

cd ..
mkdir aligned

hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa Daphnia_pulex.allmasked

#Loop through all forward and reverse paired reads and run hisat2 on each pair
for f1 in *pairedForward.fq.gz; do
	hisat2 -q -x Daphnia_pulex.allmasked -1 $f1 -2 "${f1:0:${#f1}-19}"pairedReverse.fq.gz -S aligned --aligned/"${f1:0:${#f1}-19}"/summary_file.txt	
done