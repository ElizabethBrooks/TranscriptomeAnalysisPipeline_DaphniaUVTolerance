#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveFeaturesMerged_jobOutput

# script to retrieve features from a referene fasta using a gff
# usage: bash retrieveFeaturesMerged_gffread.sh

#Required modules for ND CRC servers
module load bio/2.0

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/features_gffread"
mkdir $outFolder
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# set inputs folder
inputsPath=$inputsPath"/variantsConsensus"

# retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# retrieve genome features for the consensus
consensusFeatures=$inputsPath"/"$type"_OLYM.gff"

# retrieve file name of reference
refTag=$(basename $genomeFile)

# name output file of inputs
inputOutFile=$outFolder"/retrieveCDS_summary.txt"

# add version to output file of inputs
#gffread --version > $inputOutFile

# TO-DO
# consider adding length to the contigs of the consensus
# https://groups.google.com/g/bedtools-discuss/c/mmBumZEmd4U
# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html

# status message
echo "Generating features..."

# retrieve all cds and output translated proteins
# Pulex
gffread -C -y $outFolder"/"$refTag"_longest.pep.fa" -x $outFolder"/"$refTag"_longest.cds.fa" -g $genomeFile $genomeFeatures
# OLYM
gffread -C -y $outFolder"/"$type"_consensus_longest.pep.fa" -x $outFolder"/"$type"_consensus_longest.cds.fa" -g $inputsPath"/"$type"_consensus.fa" $consensusFeatures

# status message
echo "Features generated!"
