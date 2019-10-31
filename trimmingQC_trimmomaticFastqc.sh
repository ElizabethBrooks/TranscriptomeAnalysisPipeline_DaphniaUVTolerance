#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N trimmingQC_trimmomaticFastqc_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/trimmomatic/0.32
#Prepare for adapter trimming and quality control
cd ..
#Initialize variables
qcCountStart=0
qcCountEnd=0
score=0
dirFlag=0
runNum=0
#Retrieve input paired reads path and adapter path
readPath=$(head -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/pairedReadsPath.txt")
adapterPath=$(head -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/adapterPath.txt")
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	mkdir trimmed_run"$runNum"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for $runNum run of trimming..."
	fi
done
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*1.fq.gz; do
	#Trim extension from current file name
	curFile=$(echo $f1 | sed 's/.\.fq\.gz//')
	#Trim file path from current file name
	curFileNoPath=$(basename $f1)
	curFileNoPath=$(echo $curFileNoPath | sed 's/.\.fq\.gz//')
	#Quality control using fastqc on the first read file
	if [ "$qcCountStart" -eq 0 ]; then
		fastqc $f1 --extract
		#Determine phred score for trimming
		if grep -iF "Illumina 1.5" "$curFile"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Illumina 1.9" "$curFile"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "ERROR: Illumina encoding not found... exiting"
			exit 1
		fi
		echo "${f1:0:${#f1}-7} phred score is $score."
		#QC the first read file
		#...in progress...
		if grep -iF "WARN" "$curFile"1_fastqc/summary.txt; then
			grep -iF "WARN" "$curFile"1_fastqc/summary.txt > trimmed_run"$runNum"/"$curFileNoPath"fastqc_report.txt
		fi
		if grep -iF "FAIL" "$curFile"1_fastqc/summary.txt; then
			grep -iF "FAIL" "$curFile"1_fastqc/summary.txt > trimmed_run"$runNum"/"$curFileNoPath"fastqc_report.txt
		fi
		#Only QC one file
		qcCountStart=1
	fi
	#Perform adapter trimming on paired reads
	#using 8 threads
	trimmomatic PE -threads 8 -phred"$score" $f1 "$curFile"2.fq.gz trimmed_run"$runNum"/"$curFileNoPath"pForward.fq.gz trimmed_run"$runNum"/"$curFileNoPath"uForward.fq.gz trimmed_run"$runNum"/"$curFileNoPath"pReverse.fq.gz trimmed_run"$runNum"/"$curFileNoPath"uReverse.fq.gz ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13
	#Final quality control check using fastqc on the first trimmed paired read file
	if [ qcCountEnd = 0 ]; then
		#QC paired forward read
		#...in progress...
		fastqc trimmed_run"$runNum"/"$curFileNoPath"pForward.fq.gz --extract
		if grep -iF "WARN" trimmed_run"$runNum"/"$curFileNoPath"pForward_fastqc/summary.txt; then
			grep -iF "WARN" trimmed_run"$runNum"/"$curFileNoPath"pForward_fastqc/summary.txt > trimmed_run"$runNum"/"$curFileNoPath"fastqc_report.txt
		fi
		if grep -iF "FAIL" trimmed_run"$runNum"/"$curFileNoPath"pForward_fastqc/summary.txt; then
			grep -iF "FAIL" trimmed_run"$runNum"/"$curFileNoPath"pForward_fastqc/summary.txt > trimmed_run"$runNum"/"$curFileNoPath"fastqc_report.txt
		fi
		#Only QC one file
		qcCountEnd=1
	fi
done
