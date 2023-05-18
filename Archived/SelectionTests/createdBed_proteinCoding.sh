#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N createBed_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub createdBed_proteinCoding.sh

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

# set paths for protein coding sequence lists
geneList=$outFolder"/"$refTag"_proteinCoding_genes.txt"
transList=$outFolder"/"$refTag"_proteinCoding_transcripts.txt"

# retrieve consensus genome
consPath=$outFolder"/"$type"_consensus.fa"

# set tmp and final output files for bed12 info
totalBed=$outFolder"/"$refTag".cds.bed12"
tmpFirstBed=$outFolder"/"$refTag"_first.cds.tmp.bed12"
tmpSecondBed=$outFolder"/"$refTag"_second.cds.tmp.bed12"

# pre-clean up
rm $totalBed

# create list of protein coding sequence gene names
cat $genomeFeatures | grep -w "gene_biotype=protein_coding" | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2 > $geneList

# loop over each gene name and create a tab delimeted BED12 file
# chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
# required required required optional optional optional optional ignored ignored ignored ignored ignored
while IFS= read -r gene; do
	# status message
	echo "Processing $gene ..."
	# retrieve transcript IDs
	cat $genomeFeatures | grep -w "$gene" | awk '$3 == "mRNA"' | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2 > $transList
	# loop over each transcript for the current gene
	while IFS= read -r trans; do
		# create string with name tag and a score placeholder
		initialTags=$(echo -e "$gene\t0")
		# create string with placeholders for thickStart thickEnd itemRgb blockCount blockSizes blockStarts
		finalTags=$(echo -e "0\t0\t0\t0\t0\t0")
		# retrieve and add chrom name
		cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f1 > $tmpFirstBed
		# retrieve and add start coordinate
		cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f4 | paste $tmpFirstBed - > $tmpSecondBed
		# retrieve and add end coordinate
		cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f5 | paste $tmpSecondBed - > $tmpFirstBed
		# retrieve and add strand with initial and final tags
		cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f7 | sed "s/^/$initialTags\t/g" | sed "s/$/\t$finalTags/g" | paste $tmpFirstBed - > $tmpSecondBed
		# add info for current transcript to the total gene bed12 file
		cat $tmpSecondBed >> $totalBed
	done < $transList
done < $geneList

# copy end column to thickStart and thickEnd columns
tmpEnd=$outFolder"/"$refTag"_end.cds.tmp.bed12"
cat $totalBed | cut -f3 > $tmpEnd
cat $totalBed | cut -f1-6 | paste - $tmpEnd | paste - $tmpEnd > $tmpFirstBed
cat $totalBed | cut -f9-12 | paste $tmpFirstBed - > $tmpSecondBed
mv $tmpSecondBed $totalBed

# clean up
rm $tmpFirstBed
rm $tmpSecondBed
rm $tmpEnd
rm $geneList
rm $transList

# status message
echo "Analysis complete!"
