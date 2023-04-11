#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N consensusGenotype_jobOutput

# script to generate a consensus sequence using filtered called variants
# usage: qsub generateConsensusGenotype_bcftools.sh genotype
# usage Ex: qsub generateConsensusGenotype_bcftools.sh E05
# usage Ex: qsub generateConsensusGenotype_bcftools.sh R2
# usage Ex: qsub generateConsensusGenotype_bcftools.sh Y05
# usage Ex: qsub generateConsensusGenotype_bcftools.sh Y023

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# retrieve input genotype
genotype=$1

# set input bam file type
type="filteredMapQ"

#Make output folder
outFolder=$inputsPath"/variantsConsensus"
mkdir $outFolder

# set inputs directory name
inputsPath=$inputsPath"/variantsMerged_"$type

#Name output file of inputs
inputOutFile=$outFolder"/consensusGenotype_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

# generate Olympics consensus sequence
echo "Generating $genotype consensus..."
cat $genomeFile | bcftools consensus -s $genotype $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_"$genotype".fa"
echo "cat "$genomeFile" | bcftools consensus -s "$genotype" "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_"$genotype".fa" >> "$inputOutFile"

# index consensus fasta
samtools faidx $outFolder"/"$type"_consensus_"$genotype".fa"
