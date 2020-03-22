#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_trinity_jobOutput
#$ -pe smp 8
#Script to perform Trinity de novo transcriptome assembly
#Usage: qsub assembly_trinity.sh trimmedFolder
#Usage Ex: qsub assembly_trinity.sh trimmed_run1

#Required modules for ND CRC servers
module load bio/2.0
#Retrieve aligned reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
#Retrieve assembly outputs absolute path
outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Create output directory
outputFolder="$outputsPath"/"$1"_assembly
mkdir "$outputFolder"
#Move to outputs directory
cd "$outputFolder"
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
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
if [[ "$1"  != sortedCoordinate* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_assembly_summary.txt
echo "Sample $curSampleNoPath is being assembled..."
#Run trinity assembly with each forward and revered reads, using 8 threads
Trinity --seqType fq --max_memory 50G --left "$inputsPath"/"$1"/*_pForward.fq.gz --right "$inputsPath"/"$1"/*_pReverse.fq --CPU 8
echo "Sample $curSampleNoPath has been assembled!"
#Add run inputs to output summary file
echo "$curSampleNoPath" >> "$inputOutFile"
echo "Trinity --seqType fq --max_memory 50G --left" "$inputsPath"/"$1"/*_pForward.fq.gz "--right" "$inputsPath"/"$1"/*_pReverse.fq --CPU 8 >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"