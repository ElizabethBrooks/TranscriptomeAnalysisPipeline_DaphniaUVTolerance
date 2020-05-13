#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N uniqueMerge_fasta_jobOutput
#$ -pe smp 8
#Script to perform merge multifasta files and retain only
#the specified unique data (by sequence, ID, or both)
#Usage: qsub uniqueMerge_fasta.sh mergeBy sortedFolder genotypes
#Usage Ex: qsub uniqueMerge_fasta.sh sequence sortedCoordinate_samtoolsHisat2_run1 Y05 Y023_5 E05 R2 PA Sierra

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
	echo "ERROR: The sorted "$1" folder of fasta files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve assembly outputs absolute path
outputsPath=$(grep "multiFASTA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/multiFASTA://g")
#Retrieve selected fasta files
counter=1
fastaList=""
genotypeList=""
#Loop through all input sets of treatments and perform selected analsysis
for i in "$@"; do
	#Skip first two arguments
	if [[ $counter -eq $# ]]; then
		fastaList="$fastaList$inputsPath/$2$i_assemblyGenomeTrinity/Trinity.fasta"
		genotypeList="$genotypeList$i"
	elif [[ $counter -ge 3 ]]; then
		fastaList="$fastaList$inputsPath/$2$i_assemblyGenomeTrinity/Trinity.fasta "
		genotypeList="$genotypeList$i "
	fi
	counter=$(($counter+1))
done
echo "Fastas : $fastaList"
echo "Genotypes : $genotypeList"
#Create output directory
outputFolder="$outputsPath"
#Move to outputs directory
cd "$outputFolder"
#Name output file and folder
multiFastaFile="$outputFolder/$2$genotypeList_assemblyGenomeTrinity_multiFasta.fasta"
#Merge the set of fasta files
#echo "Beginning merging..."
#Determine which method to merge fasta files by
#if [[ "$1" == sequence ]]; then
	#Sequence identical merge
	#awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		#(FNR==1){next}
		#{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		#!(seen[seq]++){ print ">" $0 }' $fastaList > $multiFastaFile
#elif [[ "$1" == sequenceName ]]; then
	#First part of sequence name identical merge
	#awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		#(FNR==1){next}
		#{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		#{ key=substr(name,1,index(s,"|")) }
		#!(seen[key]++){ print ">" $0 }' $fastaList > $multiFastaFile
#elif [[ "$1" == sequenceAndName ]]; then
	#Sequence name and sequence identical merge
	#awk 'BEGIN{RS=">"; FS="\n"; ORS=""}
		#(FNR==1){next}
		#{ name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) }
		#!(seen[name,seq]++){ print ">" $0 }' $fastaList > $multiFastaFile
#else
	#echo "Selected merge format for fasta files not valid... exiting!"
	#exit 1
#fi