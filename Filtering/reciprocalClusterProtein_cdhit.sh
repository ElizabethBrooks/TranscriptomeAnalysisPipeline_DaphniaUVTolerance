#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N reciprocalClusterProtein_cdhit_jobOutput
#$ -pe smp 8
#Script to use blastp to translate the nucleotide sequences of a reference genome
# for searching a protein database
#Usage: qsub reciprocalClusterProtein_cdhit.sh transcriptomeFasta clusterPercent
#Usage Ex: qsub reciprocalClusterProtein_cdhit.sh trimmed_run1E05_assemblyTrinity 0.90

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "cdhitPackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/cdhitPackage://g")
#Retrieve genome reference absolute paths for querying
proteinDB=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
#Determine input database for clustering
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
else
	#Error message
	echo "Invalid protein set entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Set inputs path
inputsPath="$assemblyPath"/"$1"
#Set outputs absolute path
outputFolder="$assemblyPath"/"$1"/reciprocalClusteredProteins_cdhit_"$2"
#Set DBs of transcriptome
inputProteinDB=$(echo "$inputsPath"/decoded_transdecoder/*.transdecoder.pep)
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
inputOutFile="$outputFolder"/reciprocalClusteredProteins_"$2"_cdhit_"$1"_summary.txt
#Convert to single line fasta
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$inputProteinDB" > "$inputProteinDB".AA
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' "$proteinDB" > "$proteinDB".AA
#Use cdhit-est to search a nucelotide database
echo "Beginning cdhit-est database search..."
proteinName=$(basename "$proteinDB")
"$softsPath"/cd-hit-est-2d -i "$inputProteinDB".AA -i2 "$proteinDB".AA -o "$proteinName"_novel -c "$2" -n 10 -d 0 -M 16000 -T 8
echo "$softsPath"/cd-hit-est-2d -i "$inputProteinDB".AA -i2 "$proteinDB".AA -o "$proteinName"_novel -c "$2" -n 5 -d 0 -M 16000 -T 8 >> "$inputOutFile"
