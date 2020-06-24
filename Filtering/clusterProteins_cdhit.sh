#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N clusterProteins_cdhit_jobOutput
#$ -pe smp 8
#Script to perform clustering of sequences using cdhit
#the specified unique data (by sequence, ID, or both)
#Usage: qsub clusterProteins_cdhit.sh mergeBy clusterPercent sortedFolder genotypes
#Usage Ex: qsub clusterProteins_cdhit.sh trimmed_run1E05_assemblyTrinity 0.95
#Usage Ex: qsub clusterProteins_cdhit.sh sortedCoordinate_samtoolsHisat2_run1Y05 0.98
#Alternate usage Ex: qsub clusterProteins_cdhit.sh PA42 0.95

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "cdhitPackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/cdhitPackage://g")
#Determine input query transcriptome for blastp
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$assemblyPath"/"$1"
	#Set outputs absolute path
	outputProteinFolder="$assemblyPath"/"$1"/clustered_cdhit_"$2"
	#Set DBs of transcriptome
	inputProteinPath="$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
elif [[ "$1" == PA42 ]]; then
	#Set inputs absolut paths
	inputProteinPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
	#Set outputs absolute path
	outputProteinFolder=$(dirname "$inputProteinPath")
	outputProteinFolder="$outputProteinFolder"/clustered_cdhit_"$2"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Make protein set output directory
mkdir "$outputProteinFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputProteinFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputProteinFolder"
#Name output file of inputs
inputOutFile="$outputProteinFolder"/clustered_"$2"_cdhit_"$1"_summary.txt
#Run cd-hit to cluster proteins
"$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinPath" -p 1 -n 5 -M 16000 –d 0 -T 8
echo "$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinPath" -p 1 -n 5 -M 16000 –d 0 -T 8 >> "$inputOutFile"
