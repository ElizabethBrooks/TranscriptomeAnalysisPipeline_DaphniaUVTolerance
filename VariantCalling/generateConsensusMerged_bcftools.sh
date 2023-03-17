#!/bin/bash
#Script to generate a consensus sequence using filtered called variants
#Usage: bash generateConsensusMerged_bcftools.sh
#Usage Ex: bash generateConsensusMerged_bcftools.sh

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsPath=$inputsPath"/variantsMerged_"$type
inputsDir=$inputsPath

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Make output folder
outFolder=$inputsPath

#Name output file of inputs
inputOutFile=$outFolder"/consensus_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

#Generate consensus sequence
cat $genomeFile | bcftools consensus $inputsDir"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus.fa"
echo "cat "$genomeFile" | bcftools consensus "$inputsDir"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus.fa" >> "$inputOutFile"
