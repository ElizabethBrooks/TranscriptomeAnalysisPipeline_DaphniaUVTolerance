#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N consensus_jobOutput

# script to generate a consensus sequence using filtered called variants
# usage: qsub generateConsensusMerged_bcftools.sh genotype
# usage Ex: qsub generateConsensusMerged_bcftools.sh E05
# usage Ex: qsub generateConsensusMerged_bcftools.sh R2
# usage Ex: qsub generateConsensusMerged_bcftools.sh Y05
# usage Ex: qsub generateConsensusMerged_bcftools.sh Y023

# required modules for ND CRC servers
#module load bio

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input genotype
genotype=$1

# set input bam file type
type="filteredMapQ"

# set inputs directory name
inputsPath=$inputsPath"/variantsMerged_"$type

# set input sample list
sampleList="../InputData/fileList_Olympics_genotypes.txt"

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Make output folder
outFolder=$inputsPath

#Name output file of inputs
inputOutFile=$outFolder"/consensus_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

# generate Olympics consensus sequence
echo "Generating $genotype consensus..."
cat $genomeFile | bcftools consensus -s $genotype $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_olym.fa"
echo "cat "$genomeFile" | bcftools consensus -s "$genotype" "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_olym.fa" >> "$inputOutFile"

# clean up
rm $outFolder"/tmpList_olym.txt"
