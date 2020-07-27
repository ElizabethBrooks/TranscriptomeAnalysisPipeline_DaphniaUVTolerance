#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash searchRBH_annotations.sh annotationMethod uniqueHitsFile annotationFile
#Usage Ex: bash searchRBH_annotations.sh PANNZER /home/mae/Documents/RNASeq_Workshop_ND/reciprocalSearched_blastp/trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt /home/mae/Documents/RNASeq_Workshop_ND/AnnotationAnalysis/PA42_proteins/GO.out.txt
#Usage Ex: bash searchRBH_annotations.sh GhostKOALA /home/mae/Documents/RNASeq_Workshop_ND/reciprocalSearched_blastp/trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt /home/mae/Documents/RNASeq_Workshop_ND/AnnotationAnalysis/PA42_proteins/user_ko.txt

#Set input file paths
tail -n +2 "$2" > tmp1.txt
annotationFile="$3"

#Replace whitespace with commas for comparisons
sed -e 's/\s\+/,/g' $annotationFile > tmp2.txt

#Set outputs
outFilePath=$(dirname "$2")
outFileAnnotations=$(echo "$2" | sed 's/\.txt//g')
outFileAnnotations="$outFileAnnotations"_"$1"_matchedAnnotations.txt
outFileUnique=$(echo "$2" | sed 's/\.txt//g')
outFileUnique="$outFileUnique"_"$1"_uniqueAnnotations.txt

#Pre-clean up
echo "query,db,annotation" > $outFileAnnotations
echo "query,db,queryHit,dbHit" > $outFileUnique

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f3," tmp2.txt; then #Match
		anno=$(grep "$f3," tmp2.txt)
		echo "$f1,$f2,$anno" >> $outFileAnnotations
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> $outFileUnique
	fi
done < tmp1.txt

#Clean up
rm tmp*.txt
