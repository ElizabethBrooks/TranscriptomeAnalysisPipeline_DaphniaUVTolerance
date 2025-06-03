#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N clusterNucleotide_cdhit_jobOutput
#$ -pe smp 8
#$ -q debug
#Script to perform clustering of sequences using cdhit
#the specified unique data (by sequence, ID, or both)
#Usage: qsub clusterNucleotide_cdhit.sh mergeBy clusterPercent
#Usage Ex: qsub clusterNucleotide_cdhit.sh trimmed_run1E05_assemblyTrinity 0.98
#Usage Ex: qsub clusterNucleotide_cdhit.sh sortedCoordinate_samtoolsHisat2_run2PA_assemblyPA42_v3.0Trinity 0.98

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "cdhitPackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/cdhitPackage://g")
#Determine input database for clustering
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Set inputs path
inputsPath="$assemblyPath"/"$1"
#Set outputs absolute path
outputNucleotideFolder="$assemblyPath"/"$1"/clusteredNucleotides_cdhit_"$2"
#Set DBs of transcriptome
inputNucleotidePath=$(echo "$inputsPath"/*.fasta)
#Make nucleotide set output directory
mkdir "$outputNucleotideFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputNucleotideFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputNucleotideFolder"
#Name output file of inputs
inputOutFile="$outputNucleotideFolder"/clusteredNucleotides_"$2"_cdhit_"$1"_summary.txt
#Run cd-hit to cluster nucleotides
"$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputNucleotidePath" -p 1 -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputNucleotidePath" -p 1 -n 10 -d 0 -M 16000 -T 8 > "$inputOutFile"
