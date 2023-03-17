#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N consensus_jobOutput

# script to generate a consensus sequence using filtered called variants
# usage: qsub generateConsensus_bcftools.sh
# usage Ex: qsub generateConsensus_bcftools.sh

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
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

# retrieve list of non-tolerant genotype samples
head -2 $sampleList > $outFolder"/tmpList_nTol.txt"
# retrieve list of tolerant genotype samples
tail -2 $sampleList > $outFolder"/tmpList_tol.txt"

# generate Olympics consensus sequence
echo "Generating Olympics consensus..."
cat $genomeFile | bcftools consensus $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_olym.fa"
echo "cat "$genomeFile" | bcftools consensus "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_olym.fa" >> "$inputOutFile"

# generate non-tolerant genotypes consensus sequence
echo "Generating non-tolerant genotypes consensus..."
cat $genomeFile | bcftools consensus -S $outFolder"/tmpList_nTol.txt" $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_nTol.fa"
echo "cat "$genomeFile" | bcftools consensus -s E05,R2 "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_nTol.fa" >> "$inputOutFile"

# generate tolerant genotypes consensus sequence
echo "Generating tolerant genotypes consensus..."
cat $genomeFile | bcftools consensus -S $outFolder"/tmpList_tol.txt" $inputsPath"/"$type"_calls.flt-norm.bcf" > $outFolder"/"$type"_consensus_tol.fa"
echo "cat "$genomeFile" | bcftools consensus -s Y05,Y023 "$inputsPath"/"$type"_calls.flt-norm.bcf > "$outFolder"/"$type"_consensus_tol.fa" >> "$inputOutFile"

# clean up
rm $outFolder"/tmpList_olym.txt"
