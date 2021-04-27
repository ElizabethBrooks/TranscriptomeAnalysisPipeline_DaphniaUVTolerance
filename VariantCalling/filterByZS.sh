#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput
#Script to perform variant calling
#Usage: qsub variantCalling_samtools.sh sortedNameFolder analysisTarget
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run3 genome
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run1 trimmed_run1E05_assemblyTrinity
#Usage Ex: qsub variantCalling_samtools.sh sortedName_samtoolsHisat2_run1 sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity

#Required modules for ND CRC servers
module load bio
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
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

#Keep only unique read alignments using ZS tag
for f in "$inputsDir"/*/accepted_hits.bam; do 
	echo "Processing file $f"; path=$(dirname $f); samtools view -h -f 0x2 $f | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > "$path"/filteredZS.bam
done
