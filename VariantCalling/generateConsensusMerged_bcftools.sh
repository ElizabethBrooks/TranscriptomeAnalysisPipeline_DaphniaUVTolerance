#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N consensusMerged_jobOutput

# script to generate a consensus sequence using filtered called variants
# usage: qsub generateConsensusMerged_bcftools.sh

#Required modules for ND CRC servers
module load bio/2.0

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

#Retrieve input bam file type
type="filteredMapQ"

#Make output folder
outFolder=$inputsPath"/variantsConsensus"
mkdir $outFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# set inputs directory name
inputsPath=$inputsPath"/variantsMerged_"$type

#Name output file of inputs
inputOutFile=$outFolder"/consensusMerged_summary.txt"

#Add version to output file of inputs
bcftools --version > $inputOutFile

# generate merged consensus sequence
echo "Generating merged consensus..."
cat $genomeFile | bcftools consensus -c $outFolder"/"$type"_consensus.chain" $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus.fa"
echo "cat "$genomeFile" | bcftools consensus -c "$outFolder"/"$type"_consensus.chain "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus.fa" >> "$inputOutFile"

# index consensus fasta
samtools faidx $outFolder"/"$type"_consensus.fa"
