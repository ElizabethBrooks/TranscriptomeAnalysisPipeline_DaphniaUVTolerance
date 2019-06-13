#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2
#$ -pe smp 8

cd ..
mkdir aligned_hisat2

#Prepare for mapping
module load bio/hisat2/2.1.0
hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_hisat2/Daphnia_pulex.allmasked

#Store forward and reverse reads in arrays
FORWARDARRAY=(trimmed/*pForward.fq.gz)
REVERSEARRAY=(trimmed/*pReverse.fq.gz)

#Run hisat2 on reads stored in arrays
hisat2 -q -p 8 -x aligned_hisat2/Daphnia_pulex.allmasked -1 ${FORWARDARRAY[*]} -2 ${REVERSEARRAY[*]} -S aligned_hisat2 --summary-file aligned_hisat2/alignedSummary.txt
