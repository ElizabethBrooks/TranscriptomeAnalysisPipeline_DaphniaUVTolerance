#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveCDS_jobOutput

# script to retrieve features from a referene fasta using a gff
# usage: bash retrieveFeatures_gffread.sh

#Required modules for ND CRC servers
#module load bio

#Retrieve sorted reads input absolute path
#inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/"
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

#Retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/features_gffread"
mkdir $outFolder

# set inputs folder
inputsPath=$inputsPath"/variantsMerged_"$type

#Retrieve genome reference absolute path for alignment
#genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
genomeFile="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna"
#Retrieve genome features absolute path for alignment
#genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
genomeFeatures="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff"

# retrieve file name of reference
refTag=$(basename $genomeFile)

#Name output file of inputs
inputOutFile=$outFolder"/retrieveCDS_summary.txt"

#Add version to output file of inputs
bedtools --version > $inputOutFile

# TO-DO
# consider adding length to the contigs of the consensus
# https://groups.google.com/g/bedtools-discuss/c/mmBumZEmd4U
# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html

# status message
echo "Generating features for the reference genome..."

# retrieve all cds for the reference and output translated proteins
gffread -v -C -y $outFolder"/"$refTag"_pep.fa" -x $outFolder"/"$refTag"_cds.fa" -g $genomeFile $genomeFeatures

# retrieve all cds, discarding shorter duplicates, and output translated proteins
gffread -v -C -M -K -d $outFolder"/"$refTag"_duplicateInfo.txt" -y $outFolder"/"$refTag"_longest_pep.fa" -x $outFolder"/"$refTag"_longest_cds.fa" -g $genomeFile $genomeFeatures


# status message
echo "Generating features for the consensus genome..."

# retrieve all cds for the consensus and output translated proteins
gffread -v -C -y $outFolder"/"$type"_consensus_pep.fa" -x $outFolder"/"$type"_consensus_cds.fa" -g $inputsPath"/"$type"_consensus.fa" $genomeFeatures

# retrieve all cds, discarding shorter duplicates, and output translated proteins
gffread -v -C -M -K -d $outFolder"/"$type"_consensus_duplicateInfo.txt" -y $outFolder"/"$type"_consensus_longest_pep.fa" -x $outFolder"/"$type"_consensus_longest_cds.fa" -g $inputsPath"/"$type"_consensus.fa" $genomeFeatures
