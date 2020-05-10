#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N uniqueMerge_fasta_jobOutput
#$ -pe smp 8
#Script to perform merge multifasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: qsub uniqueMerge_fasta.sh sortedFolder genotype mergeBy
#Usage Ex: qsub uniqueMerge_fasta.sh sortedCoordinate_samtoolsHisat2_run1 Sierra sequence

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
if [[ "$1"  != sortedCoordinate* ]]; then
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
#Retrieve assembly outputs absolute path
outputsPath=$(grep "assemblingWithGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingWithGenome://g")
#Create output directory
outputFolder="$outputsPath"/"$1""$2"_assemblyGenomeTrinity
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to outputs directory
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1""$2"_assemblyGenomeTrinity_summary.txt
#Merge and re-coordinate sort the set of bam files
readFiles=$(echo "$inputsPath"/"$1"/*_"$2"_*/*.bam)
echo "Beginning merging..."
#TODO: Add input check
#First part of sequence name identical merge
awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
	(FNR==1){next}
	{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
	{ key=substr(name,1,index(s,"|")) }
	!(seen[key]++){ print ">" $0 }' file1.fasta file2.fasta file3.fasta ...
#Sequence identical merge
awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
	(FNR==1){next}
	{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
	!(seen[seq]++){ print ">" $0 }' file1.fasta file2.fasta file3.fasta ...
#Sequence name and sequence identical merge
awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
	(FNR==1){next}
	{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
	!(seen[name,seq]++){ print ">" $0 }' file1.fasta file2.fasta file3.fasta ...
