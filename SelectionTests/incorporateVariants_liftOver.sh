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
genomeTag=$(basename $genomeFile | sed 's/\.fna//g')

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# set output folder
outFolder=$inputsPath"/selectionTests"
mkdir $outFolder
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outFolder directory already exsists... please remove before proceeding."
	exit 1
fi

# move to software directory
cd $softwarePath

# status message
echo "Beginnning analysis..."

# convert BCF to VCF
inputBcf=$inputsPath"/variantsMerged/"$type"_calls.flt-norm.bcf"
outputVcf=$outFolder"/"$type"_calls.flt-norm.vcf.gz"
bcftools view $inputBcf -Oz -o $outputVcf
#bcftools view /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsMerged/filteredMapQ_calls.flt-norm.bcf -Oz -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_calls.flt-norm.vcf.gz

# index VCF
bcftools index $outputVcf -t
#bcftools index /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_calls.flt-norm.vcf.gz -t

# convert VCF to BED9+ with optional extra fields
vcfBedPrefix=$outFolder"/"$type"_OLYM"
./vcfToBed $outputVcf $vcfBedPrefix
#./UTIL/vcfToBed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_calls.flt-norm.vcf.gz /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM

## prepare gff file
##genePredFile=$outFolder"/"$genomeTag".gp"
##./ldHgGene -out=$genePredFile database table $genomeFeatures
##./UTIL/ldHgGene -out=/scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/GCF_021134715.1.gp database table /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff

# move annotations from one assembly to another
fmtChain=$inputsPath"/variantsConsensus/"$type"_consensus.chain"
vcfBed=$vcfBedPrefix".bed"
noMap=$outFolder"/"$type"_OLYM.unMapped.bed"
./liftOver -gff $genomeFeatures $fmtChain -bedPlus=9 $vcfBed $noMap
#./UTIL/liftOver -gff /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsConsensus/filteredMapQ_consensus.chain -bedPlus=9 /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM.bed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM.unMapped.bed
##./liftOver -genePred $genePredFile $fmtChain -bedPlus=9 $vcfBed $noMap
##./UTIL/liftOver -genePred /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/GCF_021134715.1.gp /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsConsensus/filteredMapQ_consensus.chain -bedPlus=9 /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM.bed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM.unMapped.bed

# retrieve chrom sizes for checking liftOver output
refSizes=$outFolder"/"$refTag".chrom.sizes"
sampleGenome=$inputsPath"/variantsConsensus/"$type"_consensus.fa"
sampleSizes=$outFolder"/"$type"_OLYM.chrom.sizes"
./$softwarePath"/"faSize -detailed -tab $genomeFile > $refSizes
./$softwarePath"/"faSize -detailed -tab $sampleGenome > $sampleSizes
#./UTIL/faSize -detailed -tab /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna > /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/GCF_021134715.1_ASM2113471v1_genomic.chrom.sizes
#./UTIL/faSize -detailed -tab /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsConsensus/filteredMapQ_consensus.fa > /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests_test/filteredMapQ_OLYM.chrom.sizes

# status message
echo "Analysis complete!"

