#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -N adapterTrimming_trimmomatic
#$ -pe smp 1
#$ -q debug

cd ..	
mkdir trimmed

#Use up to 24 (256 GB) separate processes
N=24
#Loop through all forward and reverse reads and run trimmomatic on each pair
(
for f1 in *1.fq.gz; do
	for f2 in *2.fq.gz; do
		if [ "${f1:0:${#f1}-7}" = "${f2:0:${#f2}-7}" ]; then
			((i=i%N)); ((i++==0)) && wait
			trimmomatic PE -phred64 $f1 $f2 trimmed/"${f1:0:${#f1}-7}"_pairedForward.fq.gz trimmed/"${f1:0:${#f1}-7}"_unpairedForward.fq.gz trimmed/"${f1:0:${#f1}-7}"_pairedReverse.fq.gz trimmed/"${f1:0:${#f1}-7}"_unpairedReverse.fq.gz ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13 &
		fi
	done
done
)