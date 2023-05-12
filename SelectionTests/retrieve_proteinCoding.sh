#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveProteins_jobOutput

# script to run tests for selection for each protein sequence
# usage: qsub retrieve_proteinCoding.sh

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
outFolder=$inputsPath"/features_gffread"

# retrieve input bam file type
type="filteredMapQ"

# retrieve file name of reference
refTag=$(basename $refPath)

# set paths for protein coding sequence lists
geneList=$outFolder"/"$refTag"_proteinCoding_genes.txt"
transList=$outFolder"/"$refTag"_proteinCoding_transcripts.tmp.txt"


# create list of protein coding sequence gene names
cat $genomeFeatures | grep -w "gene_biotype=protein_coding" | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2 > $geneList

# loop over each gene name and create a tab delimeted BED12 file
# chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
# required required required optional optional optional optional ignored ignored ignored ignored ignored
while IFS= read -r line; do
	# status message
	echo "Processing $line ..."
	# retrieve transcript IDs
	cat $genomeFeatures | grep -w "$line" | awk '$3 == "mRNA"' | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2 > $transList
	# loop over each transcript for the current gene
	while IFS= read -r trans; do
		# create string with a score placeholder
		scoreTag=$(echo -e "0")
		# create string with placeholders for thickStart thickEnd itemRgb blockCount blockSizes blockStarts
		endTags=$(echo -e "0\t0\t0\t0\t0\t0")
		# retrieve chrom name
		chromName=$(cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f1)
		# retrieve start coordinate
		startCoord=$(cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f4)
		# retrieve end coordinate
		endCoord=$(cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f5)
		# retrieve CDS name
		cdsName=$(cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2)
		# retrieve strand
		# prepend the CDS score and append the thickStart thickEnd itemRgb blockCount blockSizes blockStarts
		strand=$(cat $genomeFeatures | grep -w "$trans" | awk '$3 == "CDS"' | cut -f7 | sed "s/^/$scoreTag\t/g" | sed "s/$/\t$endTags/g")
	done < $transList
done < $geneList

# clean up
rm $geneList
rm $tmpRefPep
rm $tmpConPep
rm $tmpRefNuc
rm $tmpConNuc

# status message
echo "Analysis complete!"
