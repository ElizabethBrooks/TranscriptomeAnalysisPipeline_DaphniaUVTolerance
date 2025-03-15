#!/bin/bash

# script to map gene IDs to annotations

# usage: bash mapGeneToAnnotation.sh

# retrieve input counts
counts="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned.csv"

# retrieve input annotations
annotations="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/KAP4_NCBI_functional_annotation.txt"

# setup tmp file
counts_tmp="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned_tmp.csv"
tail -n+2 $counts > $counts_tmp

# set output file
outFile="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/GeneCountsAnalyzed/Formatted/cleaned_annotated.csv"
outMap="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/email/geneToAnnotation_map.csv"

# add header
header=$(head -1 $counts)
echo "annotation,"$header > $outFile
echo "gene,annotation" > $outMap

# status message
echo "Beginning mapping..."

#Loop over first set of annotations
while read line; do
	# retrieve gene ID
	geneID=$(echo $line | cut -d"," -f1 | sed "s/gene-//g")
	# status message
	echo "Processing $geneID..."
	# map annotations
	if grep -q "$geneID" $annotations; then
		# clean annotation descriptions
		d1=$(grep "$geneID" $annotations | cut -f6 | tr -d '\n' | sed "s/,/./g" | sed "s/\///g")
		d2=$(grep "$geneID" $annotations | cut -f7 | tr -d '\n' | sed "s/,/./g" | sed "s/\///g")
		# check if the description1 is uncharacterized
		if [[ $d1 == "uncharacterized protein"* ]]; then
			#status message
			#echo "Adding $d2..."
			# swap gene ID with description 1
			echo $line | sed "s/gene-$geneID/$d2,gene-$geneID/g" >> $outFile
			# add gene ID and description 2
			echo "gene-$geneID,$d2" >> $outMap
		else
			#status message
			#echo "Adding $d1..."
			# swap gene ID with description 2
			echo $line | sed "s/gene-$geneID/$d1,gene-$geneID/g" >> $outFile
			# add gene ID and description 2
			echo "gene-$geneID,$d1" >> $outMap
		fi
	else
		# add origional data
		echo "gene-$geneID,"$line >> $outFile
		# add gene ID and gene ID
		echo "gene-$geneID,gene-$geneID" >> $outMap
	fi
done < $counts_tmp

# status message
echo "Mapping complete!"

# clean up
rm $counts_tmp
