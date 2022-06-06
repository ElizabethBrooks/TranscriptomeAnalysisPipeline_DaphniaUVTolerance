#!/bin/bash
#Script to generate a consensus sequence using filtered called variants
#Usage: bash generateConsensusMerged_bcftools.sh sortedNameFolder analysisTarget
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run1 genome filteredMapQ
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredMapQ
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 genome filteredZS

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Determine what analysis method was used for the input folder of data
if [[ "$2" == *assemblyTrinity* || "$2" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == *assembly*Trinity* || "$2" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	inputsDir="$inputsPath"/"$2"/"$1"
	outputsPath="$inputsPath"/"$2"
elif [[ "$2" == genome ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	inputsDir="$inputsPath"/"$1"
	outputsPath="$inputsPath"
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi

#Make output folder
outFolder="$inputsDir"/variantCallingBcftools_"$3"
inputsDir="$inputsDir"/variantCallingBcftools_"$3"/variantsFiltered
#Name output file of inputs
inputOutFile="$outFolder"/consensus_summary.txt

#Add version to output file of inputs
bcftools --version > "$inputOutFile"

#Retrieve input bam file type
type="$3"

#Index bcf file
bcftools index "$inputsDir"/"$type"_calls.normCollapse.bcf
echo bcftools index "$inputsDir"/"$type"_calls.normCollapse.bcf >> "$inputOutFile"

#Generate consensus sequence
cat "$genomeFile" | bcftools consensus "$inputsDir"/"$type"_calls.normCollapse.bcf > "$outFolder"/"$type"_consensus.fa
echo cat "$genomeFile" | bcftools consensus "$inputsDir"/"$type"_calls.normCollapse.bcf ">" "$outFolder"/"$type"_consensus.fa >> "$inputOutFile"