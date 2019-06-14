#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -N optimizedTrimmingQC_trimmomaticFastqc
#$ -pe smp 1
#$ -N output
#$ -t 1-36:1

#Prepare for adapter trimming and quality control
#Initialize variables
qcCountStart=0
qcCountEnd=0
score=0

#Load necessary modules for ND CRC servers
module load bio
module load bio/trimmomatic/0.32

#Move to folder with .fq.gz read files
cd ..	
mkdir trimmed

#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in *1.fq.gz; do
	#Quality control using fastqc on the first read file
	if [ "$qcCountStart" -eq 0 ]; then
		fastqc $f1 --extract
		#Determine phred score for trimming
		if grep -iF "Illumina 1.5" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Illumina 1.9" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "Illumina encoding not found!"
			exit 1
		fi
		echo "${f1:0:${#f1}-7} phred score is $score."
		#QC the first read file
		#...in progress...
		if grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		if grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		#Only QC one file
		qcCountStart=1
	fi
	#Perform adapter trimming on paired reads
	trimmomatic PE -phred"$score" $f1 "${f1:0:${#f1}-7}"2.fq.gz trimmed/"${f1:0:${#f1}-7}"pForward.fq.gz trimmed/"${f1:0:${#f1}-7}"uForward.fq.gz trimmed/"${f1:0:${#f1}-7}"pReverse.fq.gz trimmed/"${f1:0:${#f1}-7}"uReverse.fq.gz ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13
	#Final quality control check using fastqc on the first trimmed paired read file
	if [ qcCountEnd = 0 ]; then
		#QC paired forward read
		#...in progress...
		fastqc trimmed/"${f1:0:${#f1}-7}"pForward.fq.gz --extract
		if grep -iF "WARN" trimmed/"${f1:0:${#f1}-7}"pForward_fastqc/summary.txt; then
			grep -iF "WARN" trimmed/"${f1:0:${#f1}-7}"pForward_fastqc/summary.txt > trimmed/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		if grep -iF "FAIL" trimmed/"${f1:0:${#f1}-7}"pForward_fastqc/summary.txt; then
			grep -iF "FAIL" trimmed/"${f1:0:${#f1}-7}"pForward_fastqc/summary.txt > trimmed/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		#Only QC one file
		qcCountEnd=1
	fi
	#Report the task number as it is completed
	echo "Task ${SGE_TASK_ID} has completed!"
done