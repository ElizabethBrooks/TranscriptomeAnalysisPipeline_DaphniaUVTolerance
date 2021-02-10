#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N clusterProtein_cdhit_jobOutput
#$ -pe smp 8
#$ -q debug
#Script to perform clustering of sequences using cdhit
#the specified unique data (by sequence, ID, or both)
#Usage: qsub clusterProtein_cdhit.sh mergeBy clusterPercent sortedFolder genotypes
#Usage Ex: qsub clusterProtein_cdhit.sh trimmed_run1E05_assemblyTrinity 0.98
#Usage Ex: qsub clusterProtein_cdhit.sh sortedCoordinate_samtoolsHisat2_run2PA_assemblyPA42_v3.0Trinity 0.98

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
outputProteinFolder="$assemblyPath"/"$1"/clusteredProteins_cdhit_"$2"
#Set DBs of transcriptome
inputProteinPath=$(echo "$inputsPath"/decoded_transdecoder/*.transdecoder.pep)
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
inputOutFile="$outputProteinFolder"/clusteredProteins_"$2"_cdhit_"$1"_summary.txt
#Convert to single line fasta
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$inputProteinPath" > "$inputProteinPath".AA
#Run cd-hit to cluster proteins
"$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinPath".AA -p 1 -n 5 -M 16000 –d 0 -T 8
echo "$softsPath"/cd-hit -o cdhit -c $2 -i "$inputProteinPath".AA -p 1 -n 5 -M 16000 –d 0 -T 8 >> "$inputOutFile"
