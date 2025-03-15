#!/bin/bash

# script to retrieve protein coding pep sequences
# usage: bash retrieve_proteinCoding.sh outFolder genomeFeatures inRefPep
# usage ex: bash retrieve_proteinCoding.sh D_pulex /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Bioinformatics/orthofinder/proteinCoding_proteomes /Users/bamflappy/PfrenderLab/ncbi_dataset_daphniaReferences/ncbi_dataset/data/D_pulex_GCF_021134715.1/genomic.gff /Users/bamflappy/PfrenderLab/ncbi_dataset_daphniaReferences/ncbi_dataset/data/D_pulex_GCF_021134715.1/D_pulex.faa


# retrieve reference tag
refTag=$1

# retrieve outputs directory
outFolder=$2

# retrieve genome features file
genomeFeatures=$3

# retrieve pep sequences file
inRefPep=$4

# set paths for protein coding sequence lists
geneList=$outFolder"/"$refTag"_proteinCoding_genes.txt"
transList=$outFolder"/"$refTag"_proteinCoding_map.txt"

# set reference multiline pep fasta to retrieve seqs
tmpRefPep=$outFolder"/"$refTag"_proteinCoding.pep.tmp.fa"
fltRefPep=$outFolder"/"$refTag"_proteinCoding.pep.fa"

# pre-clean up
rm $transList
rm $fltRefPep

# create singleline fasta of seqs
cat $inRefPep | tr -d '\n' | sed 's/>/\n>/g' > $tmpRefPep

# create list of protein coding sequence gene names
cat $genomeFeatures | grep -w "gene_biotype=protein_coding" | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2 > $geneList

# loop over each gene name
while IFS= read -r line; do
	# status message
	echo "Processing $line ..."
	# create list of protein coding sequence transcript names, and grab the first listed transcript
	transName=$(cat $genomeFeatures | grep -w "$line" | awk '$3 == "mRNA"' | head -1 | cut -f9 | cut -d ";" -f1 | cut -d "=" -f2)
	# add transcript name to the gene to transcript map file
	echo "$transName $line" >> $transList
	# reference pep fasta
	echo -en ">"$line"\t" >> $fltRefPep
	grep -w "^>$transName" $tmpRefPep | cut -d "]" -f2 >> $fltRefPep
done < $geneList

# clean up
rm $geneList
rm $tmpRefPep

# status message
echo "Analysis complete!"
