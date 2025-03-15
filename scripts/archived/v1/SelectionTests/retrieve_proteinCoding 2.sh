#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveProteins_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub retrieve_proteinCoding.sh

# load necessary modules
module load bio

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeFeatures://g")
#genomeFeatures="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")
#inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics"

# retrieve genome reference absolute path for alignment
refPath=$(grep "genomeReference" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")
#refPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# make outputs directory name
outFolder=$inputsPath"/variantsConsensus"

# retrieve input bam file type
type="filteredMapQ"

# retrieve file name of reference
refTag=$(basename $refPath)

# retrieve consensus genome
consPath=$outFolder"/"$type"_consensus.fa"

# set files for cds fasta seqs
refNuc=$outFolder"/"$refTag".cds.fa"
conNuc=$outFolder"/"$type"_consensus.cds.fa"

# set tmp and final output files for bed12 info
totalBed=$outFolder"/"$refTag".cds.bed12"

# status message
echo "Splitting reference genome file..."

# split the reference genome file and retreive gene sequences
bedtools getfasta -fi $refPath -bed $totalBed -split -name > $refNuc 

# status message
echo "Splitting consensus genome file..."

# split the consensus genome file and retreive gene sequences
bedtools getfasta -fi $consPath -bed $totalBed -split -name > $conNuc

# status message
echo "Analysis complete!"
