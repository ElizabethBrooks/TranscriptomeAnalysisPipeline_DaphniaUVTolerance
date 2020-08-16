#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_orthoFinder_jobOutput
#$ -pe smp 8
#Script to use OrthoFinder to find orthogroups and orthologs, 
# infers rooted gene trees for all orthogroups and identifies 
# all of the gene duplication events in those gene trees
#Usage: qsub search_orthoFinder.sh proteomeFasta
#Usage Ex: qsub search_orthoFinder.sh trimmed_run1E05_assemblyTrinity
#Alternate usage Ex: qsub search_orthoFinder.sh PA42_proteins

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine input query transcriptome for orthoFinder
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
	#Set outputs absolute path
	outputFolder="$assemblyPath"/"$1"/searched_orthoFinder_"$2"
elif [[ "$1" == PA42_proteins ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
	#Set outputs absolute path
	outputPath=$(dirname "$inputsPath")
	outputFolder="$outputPath"/searched_orthoFinder_"$2"	
else
	#Error message
	echo "Invalid fasta entered (species proteome expected)... exiting!"
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
inputOutFile="$outputFolder"/searched_orthoFinder_"$1"_summary.txt
#Use OrthoFinder to find orthologs
echo "Beginning OrthoFinder search..."

echo "OrthoFinder search complete!"
#Output run commands to summary file
