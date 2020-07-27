#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash searchRBH_annotations.sh transcriptomeFasta genotype PA42File outFile
#Usage Ex: bash searchRBH_annotations.sh trimmed_run1E05_assemblyTrinity E05 PA42_proteins trimmed_run1_PA42_proteins_RBH_summary.txt

#Set input file paths
hitsFile="$1"
annotationFile="$2"

#Set outputs
outFilePath=$(dirname "$1")
outFileAnnotations=$(echo "$hitsFile" | sed 's/\.txt/\_matchedAnnotations\.txt/g')
outFileUnique=$(echo "$hitsFile" | sed 's/\.txt/\_uniqueAnnotations\.txt/g')

#Pre-clean up
echo "queryHit,dbHit,db" > $outFileRBH

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f3" $annotationFile; then #Match
		echo "$f1,$f2,$f3,$f4" >> $outFileAnnotations
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> $outFileUnique
	fi
done < $hitsFile
