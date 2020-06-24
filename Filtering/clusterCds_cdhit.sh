#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N clusterCds_cdhit_jobOutput
#$ -pe smp 8
#Script to perform clustering of cds using cdhit
#Usage: qsub clusterCds_cdhit.sh mergeBy clusterPercent sortedFolder genotypes
#Usage Ex: qsub clusterCds_cdhit.sh trimmed_run1E05_assemblyTrinity 0.95
#Usage Ex: qsub clusterCds_cdhit.sh sortedCoordinate_samtoolsHisat2_run1Y05 0.98
#Alternate usage Ex: qsub clusterCds_cdhit.sh PA42 0.95

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
	outputCdsFolder="$assemblyPath"/"$1"/clusteredCds_cdhit_"$2"
	#Set DBs of transcriptome
	inputCdsPath="$inputsPath"/decoded_transdecoder/Trinity.fasta.transdecoder.cds
elif [[ "$1" == PA42 ]]; then
	#Set inputs absolut paths
	inputCdsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
	#Set outputs absolute path
	outputCdsFolder=$(dirname "$inputCdsPath")
	outputCdsFolder="$outputCdsFolder"/clusteredCds_cdhit_"$2"
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Make Cds set output directory
mkdir "$outputCdsFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputCdsFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputCdsFolder"
#Name output file of inputs
inputOutFile="$outputCdsFolder"/clusteredCds_"$2"_cdhit_"$1"_summary.txt
#Run cd-hit to cluster Cdss
"$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputCdsPath" -p 1 -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est -o cdhitEst -c $2 -i "$inputCdsPath" -p 1 -n 10 -d 0 -M 16000 -T 8 > "$inputOutFile"
