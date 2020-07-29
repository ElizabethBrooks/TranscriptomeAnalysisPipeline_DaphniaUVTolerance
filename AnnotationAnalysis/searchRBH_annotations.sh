#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N searchRBH_annotations_jobOutput
#Script to filter reciprocal blast results for best hits
#Usage: bash searchRBH_annotations.sh annotationMethod uniqueHitsFile annotationFile
#Usage Ex: bash searchRBH_annotations.sh PANNZER trimmed_run1_PA42_proteins_blastp_PA42_proteinsConsensusRBH.txt PA42_proteins/GO.out.txt
#Usage Ex: bash searchRBH_annotations.sh PANNZER trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt PA42_proteins/GO.out.txt
#Usage Ex: bash searchRBH_annotations.sh PANNZER trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt PA42_cds/GO.out.txt
#Usage Ex: bash searchRBH_annotations.sh GhostKOALA trimmed_run1_PA42_proteins_blastp_PA42_proteinsConsensusRBH.txt PA42_proteins/user_ko.txt
#Usage Ex: bash searchRBH_annotations.sh GhostKOALA trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt PA42_proteins/user_ko.txt
#Usage Ex: bash searchRBH_annotations.sh GhostKOALA trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt PA42_cds/user_ko.txt
#Usage Ex: bash searchRBH_annotations.sh Trinotate trimmed_run1_PA42_proteins_blastp_uniqueRBH.txt PA42_cds/go_annotations.txt

#Set input file paths
searchPath=$(grep "reciprocalSearch:" ../InputData/outputPaths.txt | tr -d " " | sed "s/reciprocalSearch://g")
annotationPath=$(grep "annotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/annotations://g")
outFilePath="$searchPath"/reciprocalSearched_blastp
searchFile="$outFilePath"/"$2"
annotationFile="$annotationPath"/"$3"

#Remove header
tail -n +2 $searchFile > tmp1.txt

#Replace whitespace with commas for comparisons
sed -e 's/\s\+/,/g' $annotationFile > tmp2.txt

#Set outputs
outFileAnnotations=$(basename "$searchFile" | sed 's/\.txt//g')
outFileAnnotations="$outFilePath"/"$outFileAnnotations"_"$1"_matchedAnnotations.txt
outFileUnique=$(basename "$searchFile" | sed 's/\.txt//g')
outFileUnique="$outFilePath"/"$outFileUnique"_"$1"_uniqueAnnotations.txt

#Pre-clean up
echo "query,db,dbHit,annotation" > $outFileAnnotations
echo "query,db,queryHit,dbHit" > $outFileUnique
#Determine annotation input
if [[ "$1" == "GhostKOALA" ]]; then
	echo "query,db,dbHit,annotation" > $outFileAnnotations
else
	echo "dbHit,annotation" > $outFileAnnotations
fi

#Output status message
echo "Beginning annotation search..."

#Loop over first set of annotations
while IFS=, read -r f1 f2 f3 f4
do
	#Determine annotation for DB hit
	if grep -q "$f4," tmp2.txt; then #Match
		anno=$(grep "$f4," tmp2.txt)
		#Determine annotation input
		if [[ "$1" == "GhostKOALA" ]]; then
			echo "$f1,$f2,$anno" >> $outFileAnnotations
		else
			echo "$anno" >> $outFileAnnotations
		fi
	else #Unique
		echo "$f1,$f2,$f3,$f4" >> $outFileUnique
	fi
done < tmp1.txt

#Output status message
echo "Annotation search complete!"

#Clean up
rm tmp*.txt
