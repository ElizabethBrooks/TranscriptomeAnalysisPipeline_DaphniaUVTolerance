#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2_jobOutput
#$ -pe smp 8

#Prepare for alignment
cd ..
mkdir aligned_tophat2
mkdir aligned_tophat2/out
module load bio
genomeFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/PA42.3.0.annotation.18440.gff"
#Build reference genome
bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_tophat2/Daphnia_pulex.allmasked
#Loop through all forward and reverse paired reads and run tophat2 on each pair
# using 8 threads
for f1 in trimmed/*pForward.fq.gz; do
	echo "Sample ${f1:8:${#f1}-23} is being aligned..."
	tophat2 -p 8 -G $geneomeFile -o aligned_tophat2/out/"${f1:8:${#f1}-23}" aligned_tophat2/Daphnia_pulex.allmasked $f1 "${f1:0:${#f1}-14}"pReverse.fq.gz
	echo "Sample ${f1:8:${#f1}-23} has been aligned!"
done