#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#$ -pe smp 8

#Prepare for alignment
cd ..
mkdir stats_edgeR
module load bio/python/2.7.14
module load bio/htseq/0.11.2
COUNTER=0
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in aligned_tophat2/out/*; do
        READARRAY[i]="$f1/accepted_hits.bam, "
        let COUNTER+=1
done
#Remove the last two characters (extra comma and white space) from the last element of the read file array
READARRAY[COUNTER]="$f1/accepted_hits.bam"
echo ${READARRAY}
#Use samtools to prepare mapped reads for counting
#samtools sort -T /tmp/sample/accepted_hits_sorted.bam -o aligned/sample/accepted_hits_sorted.bam aligned/sample/accepted_hits.bam
#Use htseq-count to prepare sorted reads for stats analysis in edgeR
#htseq-count -s no -m union -t gene -i trID aligned/sample/accepted_hits_sorted.bam -i genomeFilePath > sample.counts