#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1Sierra_assembly_Trinity
#Note that the genome version input is for output file naming purposes only

#Load necessary modules for ND CRC servers
module load bio/transdecoder
#module load bio/cufflinks
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
if [[ "$1"  != trimmed*assembly_Trinity ]]; then
	echo "ERROR: The "$1" folder of trimmed fq.gz files were not found... exiting"
	exit 1
fi
#Retrieve genome reference and features paths
multiFASTAPath=$(grep "assembling:" ../InputData/inputPaths.txt | tr -d " " | sed "s/assembling://g")
multiFASTA="$multiFASTAPath"/"$1"/Trinity.fasta
#Retrieve outputs absolute path
outputsPath=$(grep "decoding:" ../InputData/outputPaths.txt | tr -d " " | sed "s/decoding://g")
outFolder="$outputsPath"/decoded_"$1"
mkdir "$outFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outFolder"
#Generate your best candidate open rading frame (ORF) predictions
TransDecoder.LongOrfs -t "$multiFASTA"
#Identify peptides with homology to known proteins
#TransDecoder.Predict -t "$multiFASTA"
