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
inputsDir=$inputsPath"/variantsMerged"

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
g2gtools logo 

# set inputs
# reference genome in fasta format
REF_SEQ=$genomeFile
# sample ID
SAMPLE_ID=$type
# vcf file
VCF=$inputsPath"/variantsMerged/"$type"_calls.flt-norm.bcf"
# gene annotation file in gtf format
REF_GTF=$genomeFeatures
# number of cores
NUM_CORES=4

echo "[g2gtools::vcf2vci] Compiling variant calls..."
g2gtools vcf2vci -o ${SAMPLE_ID}.vci -s ${SAMPLE_ID} --diploid -p ${NUM_CORES} ${VCF}
#g2gtools vcf2vci -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ.vci -s filteredMapQ --diploid /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsMerged/filteredMapQ_calls.flt-norm.bcf

echo "[g2gtools::patch] Incorporating snps..."
g2gtools patch -i ${REF_SEQ} -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.patched.fa -p ${NUM_CORES}

echo "[g2gtools::transform] Incorporating indels..."
g2gtools transform -i ${SAMPLE_ID}.patched.fa -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.fa -p ${NUM_CORES}

echo "[g2gtools::convert] Lifting over gtf file..."
g2gtools convert -i ${REF_GTF} -c ${SAMPLE_ID}.vci.gz -o ${SAMPLE_ID}.gtf

echo "[g2gtools::gtf2db] Converting gtf to DB..."
g2gtools gtf2db -i ${SAMPLE_ID}.gtf -o ${SAMPLE_ID}.gtf.db

echo "[g2gtools::extract] Extracting genes..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --genes > ${SAMPLE_ID}.genes.fa

echo "[g2gtools::extract] Extracting transcripts..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --transcripts > ${SAMPLE_ID}.transcripts.fa

echo "[g2gtools::extract] Extracting exons..."
g2gtools extract -i ${SAMPLE_ID}.fa -db ${SAMPLE_ID}.gtf.db --exons > ${SAMPLE_ID}.exons.fa

echo "[g2gtools] Done." 

# close software
source deactivate
