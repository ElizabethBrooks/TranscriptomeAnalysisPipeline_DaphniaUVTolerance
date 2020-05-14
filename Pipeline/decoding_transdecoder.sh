#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#$ -pe smp 8
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1Sierra_assemblyTrinity
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyGenomeTrinity
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsTophat2_run1Sierra_assemblyGenomeTrinity

#Load necessary modules for ND CRC servers
module load bio/transdecoder
module load bio/blast+
module load bio/hmmer
module load bio/cufflinks
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
if [[ "$1"  != *assembly* ]]; then
	echo "ERROR: The $1 folder of assembly files were not found... exiting"
	exit 1
fi
#Retrieve input assembly path
inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Retrieve genome reference and features paths
multiFASTA=$(echo "$inputsPath"/"$1"/Trinity.fasta)
geneMap=$inputsPath/$1/"Trinity.fasta.gene_trans_map"
#Set outputs absolute path
outputsPath=$inputsPath/$1
outputFolder=$outputsPath/"decoded_transdecoder"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_"$outputFolder"_summary.txt
#Generate your best candidate open rading frame (ORF) predictions
echo "Beginning decoding..."
#Generate candidate ORFs
echo "Beginning transdecoder open reading frame predictions..."
TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
echo "Finished generating transdecoder open reading frame predictions!"
echo "Beginning transdecoder coding region selection..."
TransDecoder.Predict -t "$multiFASTA"
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "TransDecoder.Predict -t" "$multiFASTA" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"