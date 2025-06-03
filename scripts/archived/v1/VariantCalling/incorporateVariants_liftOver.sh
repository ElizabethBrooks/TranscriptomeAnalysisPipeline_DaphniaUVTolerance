#!/bin/bash

# script to incorporate variants and update the sample gff
# usage: bash incorporateVariants_liftOver.sh

# load necessary software
#module load bio/2.0

# set software directory
softwarePath="/afs/crc.nd.edu/user/e/ebrooks5/UTIL"

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input type
type="filteredMapQ"

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#genomeTag=$(basename $genomeFile | sed 's/\.fna//g')

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# set output folder
outFolder=$inputsPath"/variantsConsensus"

# move to software directory
cd $softwarePath

# status message
echo "Beginnning analysis..."

# convert BCF to VCF
inputBcf=$inputsPath"/variantsMerged/"$type"_calls.flt-norm.bcf"
outputVcf=$outFolder"/"$type"_calls.flt-norm.vcf.gz"
bcftools view $inputBcf -Oz -o $outputVcf

# index VCF
bcftools index $outputVcf -t

# convert VCF to BED9+ with optional extra fields
vcfBedPrefix=$outFolder"/"$type"_OLYM"
./vcfToBed $outputVcf $vcfBedPrefix

## prepare gff file
##genePredFile=$outFolder"/"$genomeTag".gp"
##./ldHgGene -out=$genePredFile database table $genomeFeatures

# move annotations from one assembly to another
fmtChain=$outFolder"/"$type"_consensus.chain"
vcfBed=$vcfBedPrefix".bed"
noMap=$outFolder"/"$type"_OLYM.unMapped.bed"
./liftOver -gff $genomeFeatures $fmtChain -bedPlus=9 $vcfBed $noMap
##./liftOver -genePred $genePredFile $fmtChain -bedPlus=9 $vcfBed $noMap

# change updated bed file extension to gff
mv $vcfBed $vcfBedPrefix".gff"

# retrieve chrom sizes for checking liftOver output
refSizes=$outFolder"/"$refTag".chrom.sizes"
sampleGenome=$outFolder"/"$type"_consensus.fa"
sampleSizes=$outFolder"/"$type"_OLYM.chrom.sizes"
./faSize -detailed -tab $genomeFile > $refSizes
./faSize -detailed -tab $sampleGenome > $sampleSizes

# status message
echo "Analysis complete!"

