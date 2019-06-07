#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2
#$ -pe smp 8

cd ..
mkdir aligned

#Prepare for mapping
module load bio/hisat2/2.1.0
hisat2-build -f /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned/Daphnia_pulex.allmasked

#Initialize counter
COUNTER=0
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in trimmed/*pairedForward.fq.gz; do
        FORWARDARRAY[i]="$f1, "
        REVERSEARRAY[i]="${f1:0:${#f1}-19}pairedReverse.fq.gz, "
        let COUNTER+=1
done

#Run hisat2 on forward and reverse reads stored in arrays
hisat2 -p 8 -q -x aligned/Daphnia_pulex.allmasked -1 ${FORWARDARRAY[*]} -2 ${REVERSEARRAY[*]} -S aligned --aligned/summary_file.txt