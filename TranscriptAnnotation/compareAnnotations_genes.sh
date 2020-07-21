#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations_genes.sh annotationFile1 annotationMatchingFile
#Usage ex: bash compareAnnotations_genes.sh PA42_cds_PA42_transcripts/PA42_cds_KO_unique.csv PA42_cds_PA42_proteins_PA42_transcripts/KO_matching.csv

echo "Beginning annotated genes comparisons..."

#Retrieve file inputs
name1=$(basename $1 | sed 's/\_.O\.csv//g')
name2=$(basename $2 | sed 's/\_.O\.csv//g')

#Name output directory
# and make, if necessary
outputPath=$(dirname $2); outputPath="$outputPath"/"genes"

#Make output directories
mkdir "$outputPath"

#Name output files
outputFile="$outputPath"/"$name1"_genes_comparisonSummary.txt
alt1File="$outputPath"/"$name1"_genes_alternate.csv
new1File="$outputPath"/"$name1"_genes_new.csv

#Add header to output files
echo "gene" > $alt1File
echo "gene" > $new1File

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1" $2; then #Alternate
		echo "$f1" >> $alt1File
	else #New
		echo "$f1" >> $new1File
	fi
done < $1

#Determine the number of annotations in each subset
total1Num=$(wc -l $1 | cut -d " " -f 1)
matchNum=$(wc -l $2 | cut -d " " -f 1)
alt1Num=$(($(wc -l $alt1File | cut -d " " -f 1)-1))
new1Num=$(($(wc -l $new1File | cut -d " " -f 1)-1))

#Output the number of annotations
echo "File Matching Uniqe New Alternate" > $outputFile
echo "$name1 $matchNum $total1Num $new1Num $alt1Num" >> $outputFile

echo "Annotated genes comparisons complete!"