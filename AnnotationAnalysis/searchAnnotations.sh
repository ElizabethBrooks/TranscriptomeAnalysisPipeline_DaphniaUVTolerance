#!/bin/bash
#$ -r n
#$ -N searchAnnotations_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: qsub searchAnnotations.sh annotationMethod uniqueHitsFile annotationFile outAnnotations outUnique
#Usage ex: qsub searchAnnotations.sh PANNZER split.txt tmp2.txt outFileAnnotations outFileUnique

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f4," "$3"; then #Match
		anno=$(grep "$f4," "$3")
		#Determine annotation input
		if [[ "$1" == "GhostKOALA" ]]; then
			echo "$f1,$f2,$anno" >> "$4"
		else
			echo "$anno" >> "$4"
		fi
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> "$5"
	fi
done < "$2"

#Clean up
rm "$2"
