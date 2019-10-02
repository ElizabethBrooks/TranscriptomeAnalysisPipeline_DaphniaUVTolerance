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
#Retrieve input read file absolute path
inputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/trimmingInput_readPath.txt"
COUNTER=0
while IFS= read -r line; do
	for word in $line; do
		#Each line contains the tags for the replicates, genotypes, or treatments
		#with each tag for the category separated by spaces
	    if [[ COUNTER -eq 0 ]]; then
	    	readFiles=$word
	    elif [[ COUNTER -eq 1 ]]; then
	    	adapterPath="$word"
	    else
	    	echo "Incorrect number of lines in statsInputs_edgeR... exiting"
	    	exit 1
	    fi
	done	
	let COUNTER+=1
done < "$inputsFile"
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
for f1 in "$readFiles"/*1.fq.gz; do
	#Quality control using fastqc on the first read file
	if [ "$qcCountStart" -eq 0 ]; then
		fastqc $f1 --extract
		#Determine phred score for trimming
		if grep -iF "Illumina 1.5" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Illumina 1.9" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "Illumina encoding not found... exiting"
			exit 1
		fi
		echo "${f1:0:${#f1}-7} phred score is $score."
		#QC the first read file
		#...in progress...
		if grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "WARN" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:0:${#f1}-7}"_fastqc_report.txt
		fi
		if grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt; then
			grep -iF "FAIL" "${f1:0:${#f1}-7}"1_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:0:${#f1}-7}"_fastqc_report.txt
		fi
		#Only QC one file
		qcCountStart=1
	fi
	#Perform adapter trimming on paired reads
	#using 8 threads
	trimmomatic PE -threads 8 -phred"$score" $f1 "${f1:0:${#f1}-7}"2.fq.gz trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pForward.fq.gz trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"uForward.fq.gz trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pReverse.fq.gz trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"uReverse.fq.gz ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13
	#Final quality control check using fastqc on the first trimmed paired read file
	if [ qcCountEnd = 0 ]; then
		#QC paired forward read
		#...in progress...
		fastqc trimmed_run"$runNum"/"${f1:0:${#f1}-7}"pForward.fq.gz --extract
		if grep -iF "WARN" trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pForward_fastqc/summary.txt; then
			grep -iF "WARN" trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pForward_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"_fastqc_report.txt
		fi
		if grep -iF "FAIL" trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pForward_fastqc/summary.txt; then
			grep -iF "FAIL" trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"pForward_fastqc/summary.txt > trimmed_run"$runNum"/"${f1:${#readFiles}:${#f1}-7}"_fastqc_report.txt
		fi
		#Only QC one file
		qcCountEnd=1
	fi
done
