#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#$ -pe smp 8

#Prepare for alignment
cd ..
#mkdir stats_edgeR
#mkdir stats_edgeR/sorted
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Loop through all forward and reverse paired reads and store the file locations in arrays
for f1 in aligned_tophat2/out/*; do
	echo "Sample ${f1:20:${#f1}-0} is being analyzed..."
	#Run samtools to prepare mapped reads for counting
	# using 8 threads
	#samtools sort -o stats_edgeR/sorted/${f1:20:${#f1}-0}/accepted_hits.sorted.bam -T /tmp/${f1:20:${#f1}-0}/accepted_hits.sorted.bam -@ 8 $f1/accepted_hits.bam
	#Run htseq-count to prepare sorted reads for stats analysis in edgeR
	#htseq-count -s no -m union -t gene -i trID $f1/accepted_hits.sorted.bam -i /afs/crc.nd.edu/group/hoth/echo_base/genome/dpulex-genepredict-v11.gff > stats_edgeR/${f1:20:${#f1}-0}.counts
done