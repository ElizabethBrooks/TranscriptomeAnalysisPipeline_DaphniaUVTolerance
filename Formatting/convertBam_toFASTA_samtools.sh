#!/bin/bash
#Script to convert bam files to fasta format
#Usage: bash convertBam_toFASTA_samtools.sh sortedSequencesFolder
#Usage Ex: bash convertBam_toFASTA_samtools.sh sortedCoordinate_samtoolsTophat2_run1

#Required modules for ND CRC servers
module load bio
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
if [[ "$1"  != sorted* ]]; then
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
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve variant calling outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
#Make directory for output fasta files
outputFolder="$outputsPath"/"$1"_converted
mkdir "$outputFolder"
#Move to outputs directory
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1"_converted_summary.txt
#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/"$1"/*/*.bam; do
	#Name of sorted and aligned file
	curAlignedSample="$f1"
	#Trim file paths from current sample folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename $curSampleNoPath)
	#Convert bam files to fasta format using samtools
	echo "Sample $curSampleNoPath fasta is being generated..."
	samtools bam2fq "$f1" | seqtk seq -A > "$outputFolder"/"$curSampleNoPath".fa
	echo "Sample $curSampleNoPath fasta has been generated!"
	#Add run inputs to output summary file
	echo "$curSampleNoPath" >> "$inputOutFile"
	echo "samtools bam2fq" "$f1" "| seqtk seq -A >" "$outputFolder""/""$curSampleNoPath"".fa" >> "$inputOutFile"
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"