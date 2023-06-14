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
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsPath"/variantsMerged_"$type

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# set output folder
outFolder=$inputsPath"/selectionTests"

# move to outputs directory
cd $outFolder

# load the required software
# g2gtools v0.1.31
source activate g2gtools

# set inputs
# reference genome in fasta format
REF=$genomeFile
# strain name (usually a column name in the vcf file), e.g., CAST_EiJ
STRAIN=$type
# vcf file for indels
VCF_INDELS=$outFolder"/"$type"_calls.flt-indel.bcf"
# vcf file for snps
VCF_SNPS=$outFolder"/"$type"_calls.flt-snp.bcf"
# gene annotation file in gtf format
GTF=$genomeFeatures

# create a chain file for mapping bases between two genomes
g2gtools vcf2chain -f ${REF} -i ${VCF_INDELS} -s ${STRAIN} -o ${STRAIN}/REF-to-${STRAIN}.chain

# patch snps onto reference genome
g2gtools patch -i ${REF} -s ${STRAIN} -v ${VCF_SNPS} -o ${STRAIN}/${STRAIN}.patched.fa

# incorporate indels onto the snp-patched genome
g2gtools transform -i ${STRAIN}/${STRAIN}.patched.fa -c ${STRAIN}/REF-to-${STRAIN}.chain -o ${STRAIN}/${STRAIN}.fa

# create custom gene annotation with respect to the new custom genome
g2gtools convert -c ${STRAIN}/REF-to-${STRAIN}.chain -i ${GTF} -f gtf -o ${STRAIN}/${STRAIN}.gtf

# create custom annotation database
g2gtools gtf2db -i ${STRAIN}/${STRAIN}.gtf -o ${STRAIN}/${STRAIN}.gtf.db

# extract regions of interest from the custom genome
g2gtools extract --transcripts -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.transcripts.fa
g2gtools extract --genes -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.genes.fa
g2gtools extract --exons -i ${STRAIN}/${STRAIN}.fa -db ${STRAIN}/${STRAIN}.gtf.db > ${STRAIN}/${STRAIN}.exons.fa

# clean up
#rm ${STRAIN}/${STRAIN}.patched.fa
#rm ${STRAIN}/${STRAIN}.patched.fa.fai
