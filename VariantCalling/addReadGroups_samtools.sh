#!/bin/bash

# script to add read groups to merged and coordinate sorted bam files
# usage: bash addReadGroups_samtools.sh sortedFolderName
# usage Ex: bash addReadGroups_samtools.sh sortedCoordinate_samtoolsHisat2_run1_merged

# required modules for ND CRC servers
#module load bio

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsDir=$inputsPath"/"$1
outputsPath=$inputsPath

# name output file of inputs
inputOutFile=$inputsDir"/addReadGroups_summary.txt"

# add version to output file of inputs
samtools --version > $inputOutFile

# loop over each genotype
for f in $inputsDir"/"*"/accepted_hits.bam"; do 
	# remove file extension
	fileOut=$(echo $f | sed 's/\.bam/_RG\.bam/g')
	# retrieve genotype
	genotype=$(dirname $f | basename)
	# add read groups using picard tools
	samtools addreplacerg -r SM:$genotype -o $fileOut $f
	# add run inputs to output summary
	echo "samtools addreplacerg -r SM:"$genotype" -o "$fileOut" "$f >> $inputOutFile
done
