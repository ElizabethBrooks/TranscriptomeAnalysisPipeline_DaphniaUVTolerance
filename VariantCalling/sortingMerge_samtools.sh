#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sortingMerge_jobOutput
#$ -pe smp 4

# script to merge alignments for each input genotype
# usage: qsub sortingMerge_samtools.sh alignedFolder genotype
# usage Ex: qsub sortingMerge_samtools.sh sortedCoordinate_samtoolsHisat2_run1 E05
# usage Ex: qsub sortingMerge_samtools.sh sortedCoordinate_samtoolsHisat2_run1 R2
# usage Ex: qsub sortingMerge_samtools.sh sortedCoordinate_samtoolsHisat2_run1 Y05
# usage Ex: qsub sortingMerge_samtools.sh sortedCoordinate_samtoolsHisat2_run1 Y023

#Required modules for ND CRC servers
#module load bio

# retrieve input folder of trimmed data
inputFolder="$1"

# retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

# retrieve input genotype
genotype="$2"

# set sorting outputs absolute path
inputsDir=$inputsPath"/"$inputFolder

# set and create merged outputs directory
outputPath=$inputsPath"/variantsCalled_samtoolsBcftools"
mkdir $outputPath

# set and create genotype outputs directory
outputFolder=$outputPath"/"$genotype
mkdir $outputFolder
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# name output file of inputs
inputOutFile="$outputPath"/"samtools_merge_summary.txt"

# add software version to output summary file
#samtools --version > $inputOutFile

# loop through all reads and sort sam/bam files for input to samtools
sampleList=$(for f1 in $inputsPath"/"*"_"$genotype"_"*"/accepted_hits.bam"; do echo $f1; done)

# merge the list of aligned samples
samtools merge -@ 4 $outputFolder"/accepted_hits.bam" $sampleList

# add inputs to output summary
echo "samtools merge -@ 4 "$outputFolder"/"accepted_hits.bam" "$sampleList >> $inputOutFile
