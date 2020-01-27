#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N qc_fastqc_jobOutput
#$ -pe smp 8
#Script to perform fastqc quality control of paired end reads
#Usage: qsub qc_fastqc.sh readsDirectory
#Usage Ex: qsub qc_fastqc.sh assembled_trinity_run1

#Required modules for ND CRC servers
module load bio
#Prepare for adapter trimming and quality control
dirFlag=0
runNum=1
#Retrieve paired reads absolute path for alignment
#TO DO: Double check read folder format for input
readPath="$1"
#Retrieve trimming outputs absolute path
outputsPath=$(grep "qc:" ../InputData/outputPaths.txt | tr -d " " | sed "s/qc://g")
#Move to outputs directory
cd "$outputsPath"
#Make a new directory for each trimming run
while [ $dirFlag -eq 0 ]; do
	qcOut=qc_fastqc_run"$runNum"
	mkdir $qcOut
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of qc..."
	fi
done
#Loop through all forward and reverse reads and run fastqc on each pair
#TO DO: Double check read folder format for input
for f1 in "$readPath"/*1.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.\.fq\.gz//')
	#Quality control check using fastqc on the first qc paired read file
	# and grep to search for warnings and fails
	if [ qcCountEnd = 0 ]; then
		#QC paired forward read
		#...in progress...
		fastqc "$qcOut"/"$curSampleNoPath"pForward.fq.gz --extract
		if grep -iF "WARN" "$qcOut"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
			grep -iF "WARN" "$qcOut"/"$curSampleNoPath"pForward_fastqc/summary.txt > "$qcOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		if grep -iF "FAIL" "$qcOut"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
			grep -iF "FAIL" "$qcOut"/"$curSampleNoPath"pForward_fastqc/summary.txt > "$qcOut"/"$curSampleNoPath"fastqc_report.txt
		fi
		#Only QC one file
		qcCountEnd=1
	fi
done
