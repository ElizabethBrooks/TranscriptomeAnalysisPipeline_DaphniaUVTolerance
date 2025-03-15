#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchTranscripts_jobOutput
#Script to use an OrthoFinder script to find the longest transcript
# variant per gene
#Usage: bash searchTranscripts_orthoFinder.sh proteomeFasta
#Usage Ex: bash searchTranscripts_orthoFinder.sh trimmed_run1E05_assemblyTrinity
#Usage Ex: bash searchTranscripts_orthoFinder.sh trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98
#Usage Ex: bash searchTranscripts_orthoFinder.sh PA42_v4.1_proteins

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Set software path
softwarePath=$(grep "orthoFinder:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/orthoFinder://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assemblyTrinity* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searchedTranscripts_orthoFinder
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$inputsPath"/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$inputsPath"/cdhit.transdecoder.pep
	else 
		inputsPath="$inputsPath"/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *assemblyGenome* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searchedTranscripts_orthoFinder
	#Determine input file type
	if [[ "$1" == */clusteredNucleotide* ]]; then
		inputsPath="$inputsPath"/cdhitEst.transdecoder.pep
	elif [[ "$1" == */clusteredProtein* ]]; then
		inputsPath="$inputsPath"/cdhit.transdecoder.pep
	else 
		inputsPath="$inputsPath"/Trinity.fasta.transdecoder.pep
	fi
elif [[ "$1" == *proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searchedTranscripts_orthoFinder
elif [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
	#Set outputs absolute path
	outputFolder="$outputPath"/searchedTranscripts_orthoFinder	
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
	outputPath=$(dirname "$inputsPath")
	inputsPath=$(echo "$outputPath"/decoded_transdecoder/*transdecoder.pep)
	#Set outputs absolute path
	outputFolder="$outputPath"/searchedTranscripts_orthoFinder	
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
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
#Name output file of inputs
inputOutFile="$outputFolder"/searchedTranscripts_orthoFinder_summary.txt
#Use OrthoFinder to find orthologs
echo "Beginning transcript search..."
#Run script to keep the longest transcript variant per gene
python "$softwarePath"/tools/primary_transcript.py "$inputsPath"
#Output run commands to summary file
echo "python "$softwarePath"/tools/primary_transcript.py "$inputsPath > "$inputOutFile"
echo "Transcript search complete!"
