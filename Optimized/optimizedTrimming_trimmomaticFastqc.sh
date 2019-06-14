#!/bin/bash
#$ -pe smp 1
#$ -N optimizedTrimmingQC_trimmomaticFastqc_jobOutput
#$ -t 1-36:1

#Prepare for adapter trimming and quality control
#Initialize variables
score=0
fileIndex=1
taskFileIndex=0
#Set the file index to be read by the current task
taskFileIndex=$((2*SGE_TASK_ID-1))
echo "Task ${SGE_TASK_ID} will trim starting at file index $taskFileIndex."
#Load necessary modules for ND CRC servers
module load bio
module load bio/trimmomatic/0.32
#Move to folder with .fq.gz read files
cd ..	
cd ..
mkdir optimizedTrimmed
#Loop through all forward and reverse reads and run trimmomatic on each pair
for f1 in *1.fq.gz; do
	#Find phred score from first task
	if [ "${SGE_TASK_ID}" -eq 1 ] && [ "$taskFileIndex" -eq "$fileIndex" ]; then
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
		#Perform adapter trimming on paired reads
		trimmomatic PE -phred"$score" $f1 "${f1:0:${#f1}-7}"2.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"pForward.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"uForward.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"pReverse.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"uReverse.fq.gz ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13
		#Report the task number as trimming is completed
		echo "Task ${SGE_TASK_ID} has completed trimming!"
		#Exit script after a single trimmomatic run
		exit 0
	elif [ "${SGE_TASK_ID}" -gt 1 ] && [ "$taskFileIndex" -eq "$fileIndex" ]; then
		#Determine phred score for trimming
		if grep -iF "Illumina 1.5" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=64
		elif grep -iF "Illumina 1.9" "${f1:0:${#f1}-7}"1_fastqc/fastqc_data.txt; then
			score=33
		else
			echo "Illumina encoding not found!"
			exit 1
		fi
		#Perform adapter trimming on paired reads
		trimmomatic PE -phred"$score" $f1 "${f1:0:${#f1}-7}"2.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"pForward.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"uForward.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"pReverse.fq.gz optimizedTrimmed/"${f1:0:${#f1}-7}"uReverse.fq.gz ILLUMINACLIP:/afs/crc.nd.edu/x86_64_linux/bio/Trimmomatic/0.32/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:13
		#Report the task number as trimming is completed
		echo "Task ${SGE_TASK_ID} has completed trimming!"
		#Exit script after a single trimmomatic run
		exit 0
	fi
	#Increment file index for task ID comparison
	fileIndex=$((fileIndex++))
done
#Report that no files were found for the current task
echo "Task ${SGE_TASK_ID} did not find files for trimming!"
