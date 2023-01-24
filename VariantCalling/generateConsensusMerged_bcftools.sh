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
inputsDir="$inputsPath"/"$1"
inputsDir=$inputsDir"/variantCallingBcftools_"$type

#Retrieve input bam file type
type="$2"

#Make output folder
outFolder=$inputsDir"/variantsConsensus"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# set inputs path
inputsDir=$inputsDir"/variantCallingBcftools_"$type"/variantsFiltered"

#Name output file of inputs
inputOutFile=$outFolder"/consensus_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

#Index bcf file
bcftools index $inputsDir"/"$type"_calls.normCollapse.bcf"
echo "bcftools index "$inputsDir"/"$type"_calls.normCollapse.bcf" >> $inputOutFile

#Generate consensus sequence
cat $genomeFile | bcftools consensus $inputsDir"/"$type"_calls.normCollapse.bcf" > $outFolder"/"$type"_consensus.fa"
echo "cat "$genomeFile" | bcftools consensus "$inputsDir"/"$type"_calls.normCollapse.bcf > "$outFolder"/"$type"_consensus.fa" >> "$inputOutFile"