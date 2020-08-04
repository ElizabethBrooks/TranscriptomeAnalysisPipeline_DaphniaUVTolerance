#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchAnnotations_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: qsub searchAnnotations.sh annotationMethod uniqueHitsFile annotationFile
#Usage ex: qsub searchAnnotations.sh split.txt tmp2.txt outFileAnnotations outFileUnique

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f4," "$2"; then #Match
		anno=$(grep "$f4," "$2")
		#Determine annotation input
		if [[ "$1" == "GhostKOALA" ]]; then
			echo "$f1,$f2,$anno" >> "$3"
		else
			echo "$anno" >> "$3"
		fi
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> "$4"
	fi
done < "$1"

#Clean up
rm "$1"
