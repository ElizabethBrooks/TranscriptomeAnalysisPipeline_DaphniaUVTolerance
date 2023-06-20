#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N incorporateVariants_jobOutput

# script to incorporate variants and update the gff
# usage: qsub incorporateVariants_liftOver_test.sh

# load necessary software
module load bio/2.0

# set software directory
softwarePath="/afs/crc.nd.edu/user/e/ebrooks5/UTIL"

# set inputs absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input type
type="filteredMapQ"

# set inputs directory name
inputsDir=$inputsPath"/variantsMerged"

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# set output folder
outFolder=$inputsPath"/selectionTests_test"
mkdir $outFolder
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# move to software directory
cd $softwarePath

# status message
echo "Beginnning file format conversions..."

# convert BCF to VCF
inputBcf=$inputsPath"/variantsMerged/"$type"_calls.flt-norm.bcf"
outputVcf=$outFolder"/"$type"_calls.flt-norm.vcf.gz"
bcftools view $inputBcf -Oz -o $outputVcf

# index VCF
bcftools index $outputVcf -t

# convert VCF to BED9+ with optional extra fields
vcfBedPrefix=$outFolder"/"$type"_OLYM"
./vcfToBed $outputVcf $vcfBedPrefix

# move annotations from one assembly to another
fmtChain=$inputsPath"/variantsConsensus/"$type"_consensus.chain"
vcfBed=$vcfBedPrefix".bed"
noMap=$outFolder"/"$type"_OLYM.unMapped.txt"
./liftOver -gff $genomeFeatures $fmtChain -bedPlus=9 $vcfBed $noMap

# status message
echo "Analysis complete!"

