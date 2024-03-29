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
inputsDir=$inputsPath"/variantsMerged_"$type

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
vcfBed=$outFolder"/"$type"_OLYM"
./vcfToBed $outputVcf $vcfBed

# format bed file
vcfBed=$outFolder"/"$type"_OLYM.bed"
fmtBed=$outFolder"/"$type"_OLYM.fmt.bed"
cat $vcfBed | cut -f 1-9 > $fmtBed

# retrieve chrom sizes
refTag=$(basename $genomeFile | sed 's/\.fna//g')
refSizes=$outFolder"/"$refTag".chrom.sizes"
./faSize -detailed -tab $genomeFile > $refSizes

# convert bed format files to psl format
vcfPsl=$outFolder"/"$type"_OLYM.psl"
./bedToPsl $refSizes $fmtBed $vcfPsl

# convert psl records to chain records
vcfChain=$outFolder"/"$type"_OLYM.chain"
./pslToChain $vcfPsl $vcfChain

## convert a GFF3 CIGAR file to a PSL file
##./gff3ToPsl [options] queryChromSizes targetChomSizes inGff3 out.psl

## tranform a psl format file to a bed format file
##./pslToBed [options] psl bed

# status message
echo "Beginnning chain file formatting..."

# reformat chain file and correct indel sizes
# map.chain has the old genome as the target and the new genome as the query
tmpIndelChain=$outFolder"/"$type"_OLYM.flt-indel.chain"
tmpSnpChain=$outFolder"/"$type"_OLYM.flt-snp.chain"
fmtChain=$outFolder"/"$type"_OLYM.fmt.chain"

# subset chain file by snp and indel lines
cat $vcfChain | grep -A1 "chain 1" | grep -v -- "^--$" > $tmpSnpChain
cat $vcfChain | grep "chain 0" > $tmpIndelChain

# clean up
rm $fmtChain

# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
# loop over each indel entry in the chain file
while IFS= read -r line; do
	# retrieve query name
	qName=$(echo $line | cut -d " " -f 8)
	# status message
	echo "Formatting $qName ..."
	# retrieve remaining chain info
	tNaSiStrStaEn=$(echo $line | cut -d " " -f 3-7)
	indelSize=$(echo $line | cut -d "/" -f 2 | cut -d " " -f 1 | wc -c)
	indelSize=$(($indelSize-1))
	qStrand=$(echo $line | cut -d " " -f 10)
	id=$(echo $line | cut -d " " -f 13)
	# update chain entries
	score=1
	tStart=$(($tEnd-1))
	qSize=$indelSize
	qStart=0
	qEnd=$indelSize
	# output updated chain line
	echo "chain $score $tNaSiStrStaEn $qName $qSize $qStrand $qStart $qEnd $id" >> $fmtChain
	echo $qSize >> $fmtChain
done < $tmpIndelChain

# add snp chain info
cat $tmpSnpChain >> $fmtChain

# clean up
#rm $tmpSnpChain
#rm $tmpIndelChain

# move annotations from one assembly to another
noMap=$outFolder"/"$type"_OLYM.unMapped.txt"
./liftOver -gff $genomeFeatures $fmtChain -bedPlus=9 $vcfBed $noMap

# status message
echo "Analysis complete!"

