#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalClusterCds_cdhit_jobOutput
#$ -pe smp 8
#Script to use cdhit to cluster the cds of a transcript set
#Usage: qsub reciprocalClusterCds_cdhit.sh transcriptomeFasta clusterPercent
#Usage Ex: qsub reciprocalClusterCds_cdhit.sh trimmed_run1E05_assemblyTrinity 0.90

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "cdhitPackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/cdhitPackage://g")
#Retrieve genome reference absolute paths for querying
cdsDB=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
#Determine input database for clustering
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
else
	#Error message
	echo "Invalid transcript set entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Set inputs path
inputsPath="$assemblyPath"/"$1"
#Set outputs absolute path
outputFolder="$assemblyPath"/"$1"/reciprocalClusteredCds_cdhit_"$2"
#Set DBs of transcriptome
inputNucleotideDB=$(echo "$inputsPath"/decoded_transdecoder/*.transdecoder.cds)
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
inputOutFile="$outputFolder"/reciprocalClusteredCds_"$2"_cdhit_"$1"_summary.txt
#Use cdhit-est to search a cds database
echo "Beginning cdhit-est database search..."
cdsName=$(basename "$cdsDB")
"$softsPath"/cd-hit-est-2d -i "$inputNucleotideDB" -i2 "$cdsDB" -o "$cdsName"_novel -c "$2" -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est-2d -i "$inputNucleotideDB" -i2 "$cdsDB" -o "$cdsName"_novel -c "$2" -n 5 -d 0 -M 16000 -T 8 >> "$inputOutFile"
