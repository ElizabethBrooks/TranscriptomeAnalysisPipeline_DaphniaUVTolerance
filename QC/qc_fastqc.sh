#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N qc_fastqc_jobOutput
#$ -pe smp 8

# Script to perform fastqc quality control of paired end reads
# Usage: qsub qc_fastqc.sh dataStage
# Usage ex: qsub qc_fastqc.sh raw
# Usage ex: qsub qc_fastqc.sh trimmed_run1

# Required modules for ND CRC servers
module load bio/2.0

# Retrieve qc outputs absolute path
outputsPath=$(grep "qc:" ../InputData/outputPaths.txt | tr -d " " | sed "s/qc://g")

# Retrieve paired reads absolute path for alignment
# TO DO: Double check read folder format for input
if [[ "$1" == "raw" ]]; then # raw data
	readPath=$(grep "pairedReads:" ../InputData/inputPaths.txt | tr -d " " | sed "s/pairedReads://g")
	# Make outputs directory
	outputFolder=$outputsPath"/qc_"$1
elif [[ "$1" == "trimmed"* ]]; then # trimmed data
	readPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
	readPath=$readPath"/"$1
	# Make outputs directory
	outputFolder=$outputsPath"/qc_"$1
fi

#Make output directory
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"

#Loop through all forward and reverse reads and run fastqc on each pair
#TO DO: Double check read folder format for input
for f1 in "$readPath"/*.fq.gz; do
	#Trim extension from current file name
	#curSample=$(echo $f1 | sed 's/.\.fq\.gz//')
	#Trim file path from current file name
	#curSampleNoPath=$(basename $f1)
	#curSampleNoPath=$(echo $curSampleNoPath | sed 's/.\.fq\.gz//')
	#Quality control check using fastqc
	# and grep to search for warnings and fails
	#QC paired forward read
	#...in progress...
	fastqc -o "$outputFolder" $f1  --extract
	#if grep -iF "WARN" "$outputFolder"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
	#	grep -iF "WARN" "$outputFolder"/"$curSampleNoPath"pForward_fastqc/summary.txt > "$outputFolder"/"$curSampleNoPath"fastqc_report.txt
	#fi
	#if grep -iF "FAIL" "$outputFolder"/"$curSampleNoPath"pForward_fastqc/summary.txt; then
	#	grep -iF "FAIL" "$outputFolder"/"$curSampleNoPath"pForward_fastqc/summary.txt >> "$outputFolder"/"$curSampleNoPath"fastqc_report.txt
	#fi
done
