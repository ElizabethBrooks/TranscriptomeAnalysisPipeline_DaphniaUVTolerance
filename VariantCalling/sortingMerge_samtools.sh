#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sortingMerge_jobOutput
#$ -pe smp 4

# script to merge alignments for each genotype
# usage: qsub sortingMerge_samtools.sh alignedFolder
# usage Ex: qsub sortingMerge_samtools.sh aligned_hisat2_run1

#Required modules for ND CRC servers
module load bio

# retrieve input folder of trimmed data
inputFolder="$1"

# retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

# set sorting outputs absolute path
inputsDir=$inputsPath"/"$inputFolder

# set and create merged outputs directory
outputPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# name output file of inputs
inputOutFile="$outputPath"/"samtools_merge_summary.txt"

# add software version to output summary file
samtools --version > $inputOutFile


## E05 genotype

# set and create genotype outputs directory
outputFolder=$outputPath"/E05"
mkdir $outputFolder

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $inputsDir"/"*"_E05_"*"/accepted_hits.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits.bam" $sampleList
echo "samtools merge -@ 4 "$outputFolder"/accepted_hits.bam "$sampleList >> $inputOutFile

## R2 genotype

# set and create genotype outputs directory
outputFolder=$outputPath"/R2"
mkdir $outputFolder

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $inputsDir"/"*"_R2_"*"/accepted_hits.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits.bam" $sampleList
echo "samtools merge -@ 4 "$outputFolder"/accepted_hits.bam "$sampleList >> $inputOutFile

## Y05 genotype

# set and create genotype outputs directory
outputFolder=$outputPath"/Y05"
mkdir $outputFolder

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $inputsDir"/"*"_Y05_"*"/accepted_hits.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits.bam" $sampleList
echo "samtools merge -@ 4 "$outputFolder"/accepted_hits.bam "$sampleList >> $inputOutFile

## Y023 genotype

# set and create genotype outputs directory
outputFolder=$outputPath"/Y023"
mkdir $outputFolder

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $inputsDir"/"*"_Y023_"*"/accepted_hits.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits.bam" $sampleList
echo "samtools merge -@ 4 "$outputFolder"/accepted_hits.bam "$sampleList >> $inputOutFile

## Oympics set

# run script to add read groups
bash addReadGroups_samtools.sh

# set and create genotype outputs directory
outputFolder=$outputPath"/OLYM"
mkdir $outputFolder

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $outputPath"/"*"/accepted_hits_readGroups.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits_readGroups.bam" $sampleList
echo "samtools merge -@ 4 "$outputFolder"/accepted_hits_readGroups.bam "$sampleList >> $inputOutFile
