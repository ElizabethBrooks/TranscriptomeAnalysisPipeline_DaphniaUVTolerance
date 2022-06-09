#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N filterSZ_jobOutput
#Script to perform bam read quaity filtering
#Usage: qsub filterByZS.sh sortedNameFolder analysisTarget
#Usage Ex: qsub filterByZS.sh sortedCoordinate_samtoolsHisat2_run1 genome
#Usage Ex: qsub filterByZS.sh sortedCoordinate_samtoolsHisat2_run3 genome

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$2" == *assemblyTrinity* || "$2" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == *assembly*Trinity* || "$2" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == genome ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsDir="$inputsPath"/"$1"
	outputsPath="$inputsPath"
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi

#Name output file of inputs
inputOutFile="$inputsDir"/mapqFiltering_summary.txt

#Add version to output file of inputs
samtools --version > "$inputOutFile"

#Keep only unique read alignments using ZS tag
for f in "$inputsDir"/*/accepted_hits.bam; do 
	echo "Processing file $f"
	path=$(dirname $f)
	samtools view -h -f 0x2 $f | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > "$path"/filteredZS.bam
	echo samtools view -h -f 0x2 $f "|" awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' "|" samtools view -h -b ">" "$path"/filteredZS.bam >> "$inputOutFile"
done
