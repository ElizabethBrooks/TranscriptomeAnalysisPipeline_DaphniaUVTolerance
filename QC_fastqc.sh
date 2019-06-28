#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -N QC_fastqc_jobOutput
#$ -pe smp 1

#Prepare for adapter trimming and quality control
#Initialize variables
qcCountStart=0
score=0
#Load necessary modules for ND CRC servers
module load bio
module load bio/trimmomatic/0.32
#Move to folder with .fq.gz read files
cd ..
dirFlag=0
runNum=0
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	mkdir trimmed_run"$runNum"
	#Check if the folder already exists
	if [ $? -ne 0 ] ; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of trimming..."
	fi
done
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in *1.fq.gz; do
	#Quality control using fastqc on the first read file
	if [ qcCountStart = 0 ]; then
		fastqc $f1 --extract
		#Determine phred score for trimming
		if grep -iF "Encoding   Illumina 1.5" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Encoding Illumina 1.9" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "Illumina encoding not found"
			exit 1
		fi
		#QC the first read file
		#...in progress...
		if grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		if grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:0:${#f1}-7}"1_fastqc_report.txt
		fi
		#Only QC one file
		qcCountStart=1
	fi
done