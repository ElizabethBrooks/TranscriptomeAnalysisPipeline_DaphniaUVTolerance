#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_featureCounts_jobOutput
#Script to perform freatureCounts counting of aligned or sorted
# paired end reads
#Usage: qsub countingDriver_featureCounts.sh assemblyFolder sortedNameFolder
#Usage Ex: qsub countingDriver_featureCounts.sh trimmed_run1E05_assemblyTrinity sortedName_samtoolsTophat2_run1

#Required modules for ND CRC servers
module load R/3.5.3
#Set paths for r script
#export PATH=/afs/crc.nd.edu/user/e/ebrooks5/R/x86_64-pc-linux-gnu-library/3.5/Rsubread/libs:$PATH
export R_LIBS=/afs/crc.nd.edu/user/e/ebrooks5/R/x86_64-pc-linux-gnu-library/3.5:$R_LIBS
#OG_DIR=$( pwd )
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine input path
if [[ "$1" == *assembly* ]]; then
	#Retrieve reads input absolute path
	outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	outputsPath="$outputsPath"/"$1"/"$2"
else
	#Error message
	echo "Invalid assembly directory entered... exiting!"
	exit 1
fi
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Move to outputs directory
outputFolder="$outputsPath"/counted_featureCounts
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Name output file of inputs
inputOutFile="$outputFolder"/counted_summary.txt
#Count reads using featureCounts via an R script
Rscript counting_featureCounts.r "$outputsPath" "$genomeFile"
#Add run to summary file
echo "Rscript counting_featureCounts.r" "$outputsPath" "$genomeFile" > "$inputOutFile"
#Copy previous summaries
cp "$outputsPath"/*.txt "$outputFolder"
