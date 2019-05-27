#!/bin/bash
#$ -M username@nd.edu
#$ -m abe
#$ -r n
#$ -N jobname
#$ -pe smp 1

mkdir aligned

bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned

#Loop through all forward and reverse paired reads and run tophat2 on each pair
for f1 in trimmed/*pairedForward.fq.gz; do
	for f2 in trimmed/*pairedReverse.fq.gz; do
		if "${f1:0:${#f1}-19}" == "${f2:0:${#f2}-19}"; then
			tophat2 -G /afs/crc.nd.edu/group/hoth/echo_base/genome/dpulex-genepredict-v11.gff -o aligned/"${f1:0:${#f1}-19}" Daphnia_pulex.allmasked $f1 $f2
		fi
	done
done