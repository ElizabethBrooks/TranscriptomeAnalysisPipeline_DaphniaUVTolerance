#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_htseq_jobOutput
#Script to perform htseq-count counting of trimmed, aligned, then name sorted
# paired end reads
#Usage: qsub counting_htseq.sh sortedNameFolder
#Usage Ex: qsub counting_htseq.sh sortedName_samtoolsHisat2_run1

#Required modules for ND CRC servers
module load bio
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir="$inputsPath"/"$1"
outputsPath="$inputsPath"

#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	outputFolder="$outputsPath"/counted_htseq_run"$runNum"
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
inputOutFile="$outputFolder"/counted_summary.txt
#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f1 in "$inputsDir"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using htseq-count
	echo "Sample $curSampleNoPath is being counted..."
	#Determine which flags to use based on sorting method
	if [[ "$1"  == sortedName* ]]; then
		#Use name sorted flag (default)
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> "$inputOutFile"
		echo "htseq-count -f bam -a 60 -s no -m union -t gene -i ID" "$curAlignedSample" "$genomeFile" ">" "$outputFolder"/"$curSampleNoPath"/counts.txt >> "$inputOutFile"
	elif [[ "$1"  == sortedCoordinate* ]]; then
		#Use coordinate sorted flag
		#https://github.com/simon-anders/htseq/issues/37
		#--secondary-alignments ignore --supplementary-alignments ignore
		#Flag to output features in sam format
		#-o "$outputFolder"/"$curSampleNoPath"/counted.sam
		htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID "$curAlignedSample" "$genomeFile" > "$outputFolder"/"$curSampleNoPath"/counts.txt
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> "$inputOutFile"
		echo "htseq-count -f bam -a 60 -r pos -s no -m union -t gene -i ID" "$curAlignedSample" "$genomeFile" ">" "$outputFolder""/""$curSampleNoPath""/counts.txt" >> "$inputOutFile"
	else
		echo "ERROR: The bam file "$f1" was not found... exiting"
		exit 1
	fi
	echo "Sample $curSampleNoPath has been counted!"
done
#Copy previous summaries
cp "$inputsDir"/*.txt "$outputFolder"

# clean up
rm -r $inputsDir
