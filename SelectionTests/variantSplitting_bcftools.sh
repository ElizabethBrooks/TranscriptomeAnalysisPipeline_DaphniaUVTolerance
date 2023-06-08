#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N variantSplitting_jobOutput

# script to perform variant splitting after calling
# usage: qsub variantSplitting_bcftools.sh

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsDir"/variantsMerged_"$type

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# set output folder
outFolder=$inputsPath"/selectionTests"

# name output file of inputs
inputOutFile=$outFolder"/variantSplitting_summary.txt"

# add version to output file of inputs
bcftools --version > $inputOutFile

# include sites with SNPs 
bcftools filter --threads 4 -i 'TYPE="snp"' $inputsDir"/"$type"_calls.flt-norm.bcf" -Ob -o $outFolder"/"$type"_calls.flt-snp.bcf"
echo "bcftools filter --threads 4 -i 'TYPE=\"snp\"' "$inputsDir"/"$type"_calls.flt-norm.bcf -Ob -o "$outFolder"/"$type"_calls.flt-snp.bcf" >> $inputOutFile
echo "& including sites with SNPs: " >> $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.flt-snp.bcf" | grep -v "#" | wc -l >> $outputsFile

# include sites with indels 
bcftools filter --threads 4 -i 'TYPE="indel"' $inputsDir"/"$type"_calls.flt-norm.bcf" -Ob -o $outFolder"/"$type"_calls.flt-indel.bcf"
echo "bcftools filter --threads 4 -i 'TYPE=\"indel\"' "$inputsDir"/"$type"_calls.flt-norm.bcf -Ob -o "$outFolder"/"$type"_calls.flt-indel.bcf" >> $inputOutFile
echo "& including sites with SNPs: " >> $outputsFile
bcftools view --threads 4 $inputsDir"/"$type"_calls.flt-indel.bcf" | grep -v "#" | wc -l >> $outputsFile
