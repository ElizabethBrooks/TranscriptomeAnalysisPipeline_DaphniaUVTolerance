#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N consensusMerged_jobOutput

# script to generate a consensus sequence using filtered called variants
# usage: qsub generateConsensusMerged_bcftools.sh
# usage Ex: qsub generateConsensusMerged_bcftools.sh

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsPath=$inputsPath"/variantsMerged_"$type

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Make output folder
outFolder=$inputsPath

#Name output file of inputs
inputOutFile=$outFolder"/consensusMerged_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

# generate merged consensus sequence
echo "Generating merged consensus..."
cat $genomeFile | bcftools consensus $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_olym.fa"
echo "cat "$genomeFile" | bcftools consensus "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_olym.fa" >> "$inputOutFile"

# clean up
rm $outFolder"/tmpList_olym.txt"
