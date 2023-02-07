#!/bin/bash
#Script to generate a consensus sequence using filtered called variants
#Usage: bash generateConsensusMerged_bcftools.sh sortedFolderName filterType
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run1 filteredMapQ
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredMapQ
#Usage Ex: bash generateConsensusMerged_bcftools.sh sortedCoordinate_samtoolsHisat2_run3 filteredZS

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="$inputsPath"/"$1"

#Retrieve input bam file type
type="$2"

# set inputs directory name
inputsPath=$inputsPath"/variantCallingMerged_"$type"_"$3
inputsDir=$inputsPath"/variantsFiltered"

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Make output folder
outFolder=$inputsPath"/variantsConsensus"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile=$outFolder"/consensus_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

#Generate consensus sequence
cat $genomeFile | bcftools consensus $inputsDir"/"$type"_calls.normCollapse.bcf" > $outFolder"/"$type"_consensus.fa"
echo "cat "$genomeFile" | bcftools consensus "$inputsDir"/"$type"_calls.normCollapse.bcf > "$outFolder"/"$type"_consensus.fa" >> "$inputOutFile"
