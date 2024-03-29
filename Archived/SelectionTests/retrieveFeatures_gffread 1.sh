#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveCDS_jobOutput

# script to retrieve features from a referene fasta using a gff
# usage: bash retrieveFeatures_gffread.sh

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/features_gffread"
mkdir $outFolder

# set inputs folder
inputsPath=$inputsPath"/variantsConsensus"

# retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#genomeFile="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#genomeFeatures="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

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
gffread -C -y $outFolder"/"$type"_consensus_longest.pep.fa" -x $outFolder"/"$type"_consensus_longest.cds.fa" -g $inputsPath"/"$type"_consensus.fa" $genomeFeatures
# E05
gffread -C -y $outFolder"/"$type"_consensus_E05_longest.pep.fa" -x $outFolder"/"$type"_consensus_E05_longest.cds.fa" -g $inputsPath"/"$type"_consensus_E05.fa" $genomeFeatures
# R2
gffread -C -y $outFolder"/"$type"_consensus_R2_longest.pep.fa" -x $outFolder"/"$type"_consensus_R2_longest.cds.fa" -g $inputsPath"/"$type"_consensus_R2.fa" $genomeFeatures
# Y05
gffread -C -y $outFolder"/"$type"_consensus_Y05_longest.pep.fa" -x $outFolder"/"$type"_consensus_Y05_longest.cds.fa" -g $inputsPath"/"$type"_consensus_Y05.fa" $genomeFeatures
# Y023
gffread -C -y $outFolder"/"$type"_consensus_Y023_longest.pep.fa" -x $outFolder"/"$type"_consensus_Y023_longest.cds.fa" -g $inputsPath"/"$type"_consensus_Y023.fa" $genomeFeatures

# status message
echo "Features generated!"
