#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput
#Script to perform htseq-count counting of trimmed, aligned, then name sorted
# paired end reads
#Usage: qsub counting_htseq.sh sortedNameFolder
#Usage Ex: qsub counting_htseq.sh sortedName_samtoolsTophat2_run1
#Required modules for ND CRC servers
module load bio
module load bio/python/2.7.14
module load bio/htseq/0.11.2
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
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
#Determine if the correct analysis folder was input
if [[ "$1"  == sortedName* ]]; then
	#Set name sorted method flag
	sortType="Name"
elif [[ "$1"  == sortedCoordinate* ]]; then
	#Set coordinate sorted method flag
	sortType="Coordinate"
else
	echo "ERROR: The "$1" folder of name or coordinate sorted files were not found... exiting"
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
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#Retrieve sorted reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "counting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/counting://g")
#Move to outputs directory
cd "$outputsPath"
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	outputFolder=counted"$sortType"_htseq"$analysisMethod"_run"$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of htseq counting of "$1" data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsPath"/"$1"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using htseq-count
	echo "Sample $curSampleNoPath is being counted..."
	#Determine which flags to use based on sorting method
	if [[ "$1"  == sortedName* ]]; then
		#Use name sorted flag (default)
		#https://github.com/simon-anders/htseq/issues/37
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -s no -m union -t gene -i ID --secondary-alignments ignore --supplementary-alignments ignore "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> $inputOutFile
		echo "htseq-count -f bam -s no -m union -t gene -i ID --secondary-alignments ignore --supplementary-alignments ignore" "$curAlignedSample" "$genomeFile" ">" "$outputFolder"/"$curSampleNoPath"/counts.txt >> $inputOutFile
	elif [[ "$1"  == sortedCoordinate* ]]; then
		#Use coordinate sorted flag
		#https://github.com/simon-anders/htseq/issues/37
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -r pos -s no -m union -t gene -i ID --secondary-alignments ignore --supplementary-alignments ignore "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> $inputOutFile
		echo "htseq-count -f bam -r pos -s no -m union -t gene -i ID --secondary-alignments ignore --supplementary-alignments ignore" "$curAlignedSample" "$genomeFile" ">" "$outputFolder"/"$curSampleNoPath"/counts.txt >> $inputOutFile
	else
		echo "ERROR: The "$1" folder of name or coordinate sorted files were not found... exiting"
		exit 1
	fi
	echo "Sample $curSampleNoPath has been counted!"
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"
