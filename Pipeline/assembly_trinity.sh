#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_trinity_jobOutput
#$ -pe smp 8
#Script to perform Trinity de novo transcriptome assembly
#Usage: qsub assembly_trinity.sh trimmedFolder genotype
#Usage Ex: qsub assembly_trinity.sh trimmed_run1 Sierra

#Required modules for ND CRC servers
module load bio/2.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != trimmed* ]]; then
	echo "ERROR: The "$1" folder of trimmed fq.gz files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
inputsFolder="$inputsPath"/"$1"
samplesPath=$(echo ../InputData/samplesFile_trinity.txt)
#Retrieve assembly outputs absolute path
outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Create output directory
outputFolder="$outputsPath"/"$1""$2"_assembly_Trinity
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to inputs directory
cd "$inputsFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1""$2"_assembly_summary.txt
#Re-set reads file paths using input genotype tag
sed "s/GENEOTYPE/$2/g" "$samplesPath" > tmpSamplesFile.txt
#Run trinity assembly with each forward and revered reads, using 8 threads
#Trinity --seqType fq --max_memory 50G --samples_file tmpSamplesFile.txt --CPU 8 --output "$outputFolder"
rm tmpSamplesFile.txt
#Add run inputs to output summary file
echo "$curSampleNoPath" >> "$inputOutFile"
echo "Trinity --seqType fq --max_memory 50G --samples_file" "$samplesPath" "--CPU 8 --output" "$outputFolder" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"