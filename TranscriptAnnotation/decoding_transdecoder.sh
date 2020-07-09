#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N decoding_transdecoder_jobOutput
#$ -pe smp 1
#Script to predict coding regions from a de novo assembled transcriptome fasta file
# using Transdecoder
#Usage: qsub decoding_transdecoder.sh deNovoAssembledTranscriptomeFolder
#Usage Ex: qsub decoding_transdecoder.sh trimmed_run1Sierra_assemblyTrinity
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyGenomeTrinity
#Usage Ex: qsub decoding_transdecoder.sh sortedCoordinate_samtoolsTophat2_run1Sierra_assemblyGenomeTrinity
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_cds
#Alternate usage Ex: qsub decoding_transdecoder.sh PA42_transcripts

#Load necessary modules for ND CRC servers
module load bio/2.0
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
#Determine input query transcriptome for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve input assembly path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve genome reference and features paths
	multiFASTA=$(echo "$inputsPath"/"$1"/Trinity.fasta)
	geneMap=$inputsPath/$1/"Trinity.fasta.gene_trans_map"
	#Set outputs absolute path
	outputsPath=$inputsPath/$1
	outputFolder=$outputsPath/"decoded_transdecoder"
elif [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	#Retrieve genome reference and features paths
	multiFASTA=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	geneMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	outputFolder=$outputsPath/"decoded_transdecoder"
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	#Retrieve genome reference and features paths
	multiFASTA=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
	geneMap=$(grep "geneTransMap:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneTransMap://g")
	#Set outputs absolute path
	outputsPath=$(dirname $inputsPath)
	outputFolder=$outputsPath/"decoded_transdecoder"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Make output folder
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_decoded_transdecoder_summary.txt
#Generate your best candidate open rading frame (ORF) predictions
echo "Beginning decoding..."
#Generate candidate ORFs
echo "Beginning transdecoder open reading frame predictions..."
TransDecoder.LongOrfs -t "$multiFASTA" --gene_trans_map "$geneMap"
#Output run commands to summary file
echo "TransDecoder.LongOrfs -t" "$multiFASTA" "--gene_trans_map" "$geneMap" > "$inputOutFile"
echo "Finished generating transdecoder open reading frame predictions!"
echo "Beginning transdecoder coding region selection..."
if [[ "$1" == PA42* ]]; then
	TransDecoder.Predict -t "$multiFASTA" --no_refine_starts
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" "--no_refine_starts" >> "$inputOutFile"
else
	TransDecoder.Predict -t "$multiFASTA"
	#Output run commands to summary file
	echo "TransDecoder.Predict -t" "$multiFASTA" >> "$inputOutFile"
fi
echo "Finished transdecoder coding region selection!"
echo "Decoding complete!"
