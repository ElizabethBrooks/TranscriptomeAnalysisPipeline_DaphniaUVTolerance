#!/bin/bash
#Script to generate a multi FASTA file from aligned transcript sequences
#Usage: bash generateMultiFASTA_bash.sh sortedSequencesFolder
#Usage Ex: bash generateMultiFASTA_bash.sh variants_samtoolsTophat2_run1

#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != variants* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder of FASTA files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "variants:" ../InputData/outputPaths.txt | tr -d " " | sed "s/variants://g")
#Retrieve outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
outFolder="$outputsPath"/"$1"_multiFASTA
mkdir "$outFolder"
#Move to outputs directory
cd "$outputsPath"
#Loop through all FASTA files and combine into a single multiFASTA file
for f1 in "$inputsPath"/"$1"/*.fa; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Combine paired-end read FASTA files into a multiFASTA
	echo "Sample $curSampleNoPath fasta is being added..."
	cat "$f1" >> "$outFolder"/"$1"_multiFASTA.fa
	echo "Sample $curSampleNoPath fasta has been added!"
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"