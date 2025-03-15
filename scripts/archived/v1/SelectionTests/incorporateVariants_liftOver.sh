#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N incorporateVariants_jobOutput

# script to incorporate variants and update the gff
# usage: qsub incorporateVariants_liftOver.sh

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
inputsDir=$inputsPath"/variantsMerged_"$type

# retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeReference" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

# set output folder
outFolder=$inputsPath"/selectionTests"

# convert BCF to VCF
inputBcf=$inputsPath"/variantsMerged/"$type"_calls.flt-norm.bcf"
outputVcf=$outFolder"/"$type"_calls.flt-norm.vcf.gz"
bcftools view $inputBcf -Oz -o $outputVcf
#bcftools view /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/variantsMerged/filteredMapQ_calls.flt-norm.bcf -Oz -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-norm.vcf.gz

# index VCF
bcftools index $outputVcf -t
#bcftools index /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-norm.vcf.gz -t

# convert VCF to BED9+ with optional extra fields
vcfBed=$outFolder"/"$type"_OLYM"
./$softwarePath"/"vcfToBed $outputVcf $vcfBed
#./UTIL/vcfToBed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-norm.vcf.gz /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM

# format bed file
vcfBed=$outFolder"/"$type"_OLYM.bed"
fmtBed=$outFolder"/"$type"_OLYM.fmt.bed"
cat $vcfBed | -f 1-9 > $fmtBed
#cat /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.bed | cut -f 1-9 > /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.fmt.bed

# retrieve chrom sizes
refTag=$(basename $genomeFile | sed 's/\.fna//g')
refSizes=$outFolder"/"$refTag".chrom.sizes"
./$softwarePath"/"faSize -detailed -tab $genomeFile > $refSizes
#./UTIL/faSize -detailed -tab /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna > /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.chrom.sizes

# convert bed format files to psl format
vcfPsl=$outFolder"/"$type"_OLYM.psl"
./$softwarePath"/"bedToPsl $refSizes $fmtBed $vcfPsl
#./UTIL/bedToPsl /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.chrom.sizes /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.fmt.bed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.psl

# convert psl records to chain records
vcfChain=$outFolder"/"$type"_OLYM.chain"
./$softwarePath"/"pslToChain $vcfPsl $vcfChain
#./UTIL/pslToChain /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.psl /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.chain

## convert a GFF3 CIGAR file to a PSL file
##./$softwarePath"/"gff3ToPsl [options] queryChromSizes targetChomSizes inGff3 out.psl
##./UTIL/gff3ToPsl /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.chrom.sizes /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.chrom.sizes /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.psl

## tranform a psl format file to a bed format file
##./$softwarePath"/"pslToBed [options] psl bed
##./UTIL/pslToBed

# reformat chain file and correct indel sizes
# map.chain has the old genome as the target and the new genome as the query
tmpIndelChain=$outFolder"/"$type"_OLYM.flt-indel.chain"
tmpSnpChain=$outFolder"/"$type"_OLYM.flt-snp.chain"
fmtChain=$outFolder"/"$type"_OLYM.fmt.chain"
cat $pslChain | grep "chain 0" > $tmpIndelChain
cat $pslChain | grep "chain 1" > $tmpSnpChain

# add snp chain info
cat $tmpSnpChain > $fmtChain

# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
# loop over each indel entry in the chain file
while IFS= read -r line; do
	# retrieve chain info
	lineStart=$(echo $line | cut -d " " -f 1-5)
	tEnd=$(echo $line | cut -d " " -f 7)
	qName=$(echo $line | cut -d " " -f 8)
	indelSize=$(echo $line | cut -d "/" -f 2 | cut -d " " -f 1 | wc -c)
	indelSize=$(($indelSize-1))
	qStrand=$(echo $line | cut -d " " -f 10)
	id=$(echo $line | cut -d " " -f 13)
	# update chain entries
	tStart=$(($tEnd-1))
	qSize=$indelSize
	qStart=0
	qEnd=$indelSize
	# output updated chain line
	echo "$lineStart $tStart $tEnd $qName $qSize $qStart $qEnd $id" >> $fmtChain
	echo $qSize >> $fmtChain
done < $tmpIndelChain

# move annotations from one assembly to another
noMap=$outFolder"/"$type"_OLYM.unMapped.txt"
./$softwarePath"/"liftOver -gff $genomeFeatures $fmtChain -bedPlus=9 $vcfBed $noMap
#./UTIL/liftOver -gff /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/genomic.gff /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.chain -bedPlus=9 /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.bed /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_OLYM.unmapped.bed


####

## activate the python2 environment for local run
##source /afs/crc.nd.edu/user/e/ebrooks5/.bashrc
##conda activate /afs/crc.nd.edu/user/e/ebrooks5/.conda/envs/python2

## load the required software
## g2gtools v0.1.31
##conda activate g2gtools

## replace lowercase character with uppercase
##cat /afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/daphnia_pulex/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna | tr '[a-z]' '[A-Z]' > /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.fa

## set inputs
## set reference genome tag
##refTag="KAP4"
## reference genome in fasta format
##ref=$genomeFile
##ref="/scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.fa"
## strain name (usually a column name in the vcf file), e.g., CAST_EiJ
##strain=$type
## vcf file for indels
##vcf_indels=$outFolder"/"$type"_calls.flt-indel.bcf"
## vcf file for snps
##vcf_snps=$outFolder"/"$type"_calls.flt-snp.bcf"
## gene annotation file in gtf format
##gtf=$genomeFeatures

## convert BCF to VCF
##bcftools view /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-indel.bcf -Oz -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-indel.vcf.gz
##bcftools view /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-snp.bcf -Oz -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-snp.vcf.gz

## index VCFs
##bcftools index /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-indel.vcf.gz -t
##bcftools index /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-snp.vcf.gz -t

## create a chain file for mapping bases between two genomes
##g2gtools vcf2chain -f $ref -i $vcf_indels -s $strain -o $strain"/"$refTag"-to-"$strain".chain"
##g2gtools vcf2chain -f /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.fa -i /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-indel.vcf.gz -s "E05" -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/KAP4-to-E05.chain

## patch snps onto reference genome
##g2gtools patch -i $ref -s $strain -v $vcf_snps -o $strain"/"$strain".patched.fa"
##g2gtools patch -i /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/GCF_021134715.1_ASM2113471v1_genomic.fa -s "E05" -v /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/filteredMapQ_calls.flt-snp.vcf.gz -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/E05.patched.fa

## incorporate indels onto the snp-patched genome
##g2gtools transform -i $strain/$strain.patched.fa -c $$strain"/"$refTag"-to-"$strain".chain" -o $strain"/consensus_"$strain".fa"
##g2gtools transform -i /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/E05.patched.fa -c /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/KAP4-to-E05.chain -o /scratch365/ebrooks5/OLYM_dMelUV_analysis/KAP4_NCBI/variantsCalled_samtoolsBcftools/selectionTests/consensus_E05.fa

## create custom gene annotation with respect to the new custom genome
##g2gtools convert -c $strain"/"$refTag"-to-"$strain".chain" -i $gtf -f gtf -o $strain/$strain.gtf

## clean up
##rm $strain"/"$strain".patched.fa"
##rm $strain"/"$strain".patched.fa.fai"

## close the software environment
##conda deactive
