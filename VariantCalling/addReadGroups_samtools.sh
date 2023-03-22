#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N addRG_jobOutput

# script to add read groups to merged and coordinate sorted bam files
# usage: qsub addReadGroups_samtools.sh

# required modules for ND CRC servers
module load bio

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/variantsCalled_samtoolsBcftools"
outputsPath=$inputsPath

# name output file of inputs
inputOutFile=$inputsDir"/addReadGroups_summary.txt"

# add version to output file of inputs
samtools --version > $inputOutFile

# loop over each genotype
for f in $inputsDir"/"*"/accepted_hits.bam"; do 
	# remove file extension
	fileOut=$(echo $f | sed 's/\.bam/_readGroups\.bam/g')
	# retrieve genotype
	genotype=$(dirname $f)
	genotype=$(basename $genotype)
	# add read groups using picard tools
	samtools addreplacerg -@ 4 -r ID:"OLYM_"$genotype -r SM:$genotype -o $fileOut $f
	# add run inputs to output summary
	echo "samtools addreplacerg -@ 4 -r SM:"$genotype" -o "$fileOut" "$f >> $inputOutFile
	# clean up
	rm $f
done
