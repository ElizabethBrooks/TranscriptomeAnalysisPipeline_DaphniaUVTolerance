#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N generateTranscriptsCufflinks_jobOutput
#Script to generate a multi FASTA file for all transcripts in a GFF file
#Usage: qsub generateTranscripts_cufflinks.sh sortedNameFolder analysisTarget
#Usage Ex: qsub generateTranscripts_cufflinks.sh sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ
#Usage Ex: qsub generateTranscripts_cufflinks.sh genomeReference

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome reference and features absolute paths
genomeFeatFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

#Determine what analysis method was used for the input folder of data
if [[ "$1" == sorted* ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	type=$(echo "$2" | cut -d"_" -f2)
	inputsPath="$inputsPath"/"$1"/"$2"/"$type"_consensus.fa
elif [[ "$1" == genomeReference ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi

#Name output file of inputs
outputsPath=$(dirname "$inputsPath")
inputOutFile="$outFolder"/generateTranscipts_summary.txt

#Generate a fasta index file using samtools
#samtools faidx "$genomeFile"

#Generate a FASTA file with the DNA sequences for all transcripts in the GFF file
gffread -w "$outFolder"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile" 
echo gffread -w "$outFolder"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile" > "$inputOutFile"