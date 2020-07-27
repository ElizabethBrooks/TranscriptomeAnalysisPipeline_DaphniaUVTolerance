#!/bin/bash
#Script to filter reciprocal blast results for best hits
#Usage: bash searchRBH_annotations.sh uniqueHitsFile annotationFile
#Usage Ex: bash searchRBH_annotations.sh 

#Set input file paths
hitsFile="$1"
annotationFile="$2"

#Set outputs
outFilePath=$(dirname "$1")
outFileAnnotations=$(echo "$hitsFile" | sed 's/\.txt/\_matchedAnnotations\.txt/g')
outFileUnique=$(echo "$hitsFile" | sed 's/\.txt/\_uniqueAnnotations\.txt/g')

#Pre-clean up
echo "annotation" > $outFileAnnotations
echo "query,db,queryHit,dbHit" > $outFileUnique

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f3" $annotationFile; then #Match
		grep "$f3" $annotationFile >> $outFileAnnotations
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> $outFileUnique
	fi
done < $hitsFile
