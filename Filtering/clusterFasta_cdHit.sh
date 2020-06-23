#!/bin/bash
#Script to perform merge multifasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: bash clusterFasta_cdHit.sh mergeBy clusterPercent sortedFolder genotypes
#Usage Ex: bash clusterFasta_cdHit.sh sortedCoordinate_samtoolsHisat2_run1Y05 0.98
#Usage Ex: bash clusterFasta_cdHit.sh sortedCoordinate_samtoolsHisat2_run1Y05 0.98
#Alternate usage Ex: bash clusterFasta_cdHit.sh PA42 0.95

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
	outputFolder="$assemblyPath"/"$1"/clustered_cdhit_"$2"
elif [[ "$1" == PA42 ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
	inputsPath=$(dirname "$inputsPath")
	#Set outputs absolute path
	outputFolder="$inputsPath"/clustered_cdhit_"$2"
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
inputOutFile="$outputFolder"/clustered_"$2"_cdhit_"$1"_summary.txt
#Run cd-hit to cluster the current transcriptome
"$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputsPath"/Trinity.fasta -p 1 -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputsPath"/Trinity.fasta -p 1 -n 10 -d 0 -M 16000 -T 8 > "$inputOutFile"
#Run cd-hit to cluster the current transcriptome protein set
"$softsPath"/cd-hit -o cdhit -c $2 -i "$inputsPath"/Trinity.fasta.transdecoder.pep -p 1 -n 5 -M 16000 –d 0 -T 8
echo "$softsPath"/cd-hit -o cdhit -c $2 -i "$inputsPath"/Trinity.fasta.transdecoder.pep -p 1 -n 5 -M 16000 –d 0 -T 8 >> "$inputOutFile"
