#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_tuxedo_jobOutput
#$ -pe smp 8

#Prepare for alignment
cd ..
mkdir stats_tuxedo
module load bio/cufflinks/2.2.1
COUNTER=0
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in aligned_tophat2/out/*; do
        READARRAY[COUNTER]="$f1/accepted_hits.bam, "
        let COUNTER+=1
done
#Re set the last array element to rmove the last two characters
# (extra comma and white space) from the last element of the read file array
unset 'READARRAY[${#READARRAY[@]}-1]'
READARRAY[COUNTER]="$f1/accepted_hits.bam"
#Run cuffdiff on the aligned reads stored in the file array using 8 threads
cuffdiff -p 8 -o stats_tuxedo /afs/crc.nd.edu/group/hoth/echo_base/genome/dpulex-genepredict-v11.gff "${READARRAY[@]}"
