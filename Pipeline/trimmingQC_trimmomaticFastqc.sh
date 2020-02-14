#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N trimmingQC_trimmomaticFastqc_jobOutput
#$ -pe smp 8
#Script to perform trimmomatic trimming and fastqc quality control
# of paired end reads
#Usage: qsub trimmingQC_trimmomaticFastqc.sh
#Usage Ex: qsub trimmingQC_trimmomaticFastqc.sh

#Required modules for ND CRC servers
module load bio
module load bio/trimmomatic/0.32
#Prepare for adapter trimming and quality control
qcCountStart=0
qcCountEnd=0
score=0
dirFlag=0
runNum=1
#Retrieve paired reads absolute path for alignment
readPath=$(grep "pairedReads:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pairedReads://g")
#Retrieve adapter absolute path for alignment
adapterPath=$(grep "adapter:" ../InputData/inputPaths.txt | tr -d " " | sed "s/adapter://g")
#Retrieve trimming outputs absolute path
outputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
#Move to outputs directory
cd "$outputsPath"
#Make a new directory for each trimming run
while [ $dirFlag -eq 0 ]; do
	trimOut=trimmed_run"$runNum"
	mkdir $trimOut
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of trimming..."
	fi
done
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in "$readPath"/*1.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.\.fq\.gz//')
	#Quality control using fastqc on the first read file
	if [ "$qcCountStart" -eq 0 ]; then
		fastqc $f1 --extract
		#Determine phred score for trimming
		if grep -iF "Illumina 1.5" "$curSample"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Illumina 1.9" "$curSample"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "ERROR: Illumina encoding not found... exiting"
			exit 1
		fi
		echo "${f1:0:${#f1}-7} phred score is $score."
		#QC the first read file
		#...in progress...
		if grep -iF "WARN" "$curSample"1_fastqc/summary.txt; then
			grep -iF "WARN" "$curSample"1_fastqc/summary.txt > "$trimOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		if grep -iF "FAIL" "$curSample"1_fastqc/summary.txt; then
			grep -iF "FAIL" "$curSample"1_fastqc/summary.txt > "$trimOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		#Only QC one file
		qcCountStart=1
	fi
	#Name output file of inputs
	inputOutFile="$trimOut"/"$trimOut"_summary.txt
	#Perform adapter trimming on paired reads
	#using 8 threads
	#removed: HEADCROP:13
	trimmomatic PE -threads 8 -phred"$score" $f1 "$curSample"2.fq.gz "$trimOut"/"$curSampleNoPath"pForward.fq.gz "$trimOut"/"$curSampleNoPath"uForward.fq.gz "$trimOut"/"$curSampleNoPath"pReverse.fq.gz "$trimOut"/"$curSampleNoPath"uReverse.fq.gz ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	#Add run inputs to output summary file
	echo $curSampleNoPath >> $inputOutFile
	echo trimmomatic PE -threads 8 -phred"$score" $f1 "$curSample"2.fq.gz "$trimOut"/"$curSampleNoPath"pForward.fq.gz "$trimOut"/"$curSampleNoPath"uForward.fq.gz "$trimOut"/"$curSampleNoPath"pReverse.fq.gz "$trimOut"/"$curSampleNoPath"uReverse.fq.gz ILLUMINACLIP:"$adapterPath" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> $inputOutFile
	#Final quality control check using fastqc on the first trimmed paired read file
	if [ qcCountEnd = 0 ]; then
		#QC paired forward read
		#...in progress...
		fastqc "$trimOut"/"$curSampleNoPath"pForward.fq.gz --extract
		if grep -iF "WARN" "$trimOut"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
			grep -iF "WARN" "$trimOut"/"$curSampleNoPath"pForward_fastqc/summary.txt > "$trimOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		if grep -iF "FAIL" "$trimOut"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
			grep -iF "FAIL" "$trimOut"/"$curSampleNoPath"pForward_fastqc/summary.txt > "$trimOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		#Only QC one file
		qcCountEnd=1
	fi
done
