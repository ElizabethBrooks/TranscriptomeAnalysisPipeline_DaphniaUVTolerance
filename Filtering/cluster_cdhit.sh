#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N clustering_cdhit_jobOutput
#$ -pe smp 8
#Script to perform clustering of sequences using cdhit
#the specified unique data (by sequence, ID, or both)
#Usage: qsub cluster_cdhit.sh mergeBy clusterPercent sortedFolder genotypes
#Usage Ex: qsub cluster_cdhit.sh trimmed_run1E05_assemblyTrinity 0.95
#Usage Ex: qsub cluster_cdhit.sh sortedCoordinate_samtoolsHisat2_run1Y05 0.98
#Alternate usage Ex: qsub cluster_cdhit.sh PA42 0.95

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
	#Set DBs of transcriptome
	inputNucleotideDB="$inputsPath"/Trinity.fasta
	inputProteinDB="$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
elif [[ "$1" == PA42 ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
	inputsPath=$(dirname "$inputsPath")
	#Set inputs absolut paths
	inputNucleotideDB=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	inputProteinDB=$(grep "transcriptomeDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptomeDB://g")
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
#Run cd-hit to cluster proteins
"$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinDB" -p 1 -n 5 -M 16000 –d 0 -T 8
echo "$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinDB" -p 1 -n 5 -M 16000 –d 0 -T 8 >> "$inputOutFile"
#Run cd-hit to cluster nucleotides
"$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputNucleotideDB" -p 1 -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputNucleotideDB" -p 1 -n 10 -d 0 -M 16000 -T 8 > "$inputOutFile"
