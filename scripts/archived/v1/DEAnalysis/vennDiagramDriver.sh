#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash vennDiagramDriver.sh analysisResultsFile
#Usage Ex: bash vennDiagramDriver.sh exactTest_topTags.csv genotypeList
#Usage Ex: bash vennDiagramDriver.sh glmQLF_2WayANOVA_topTags_filtered.csv

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve analysis inputs path
inputPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

#Determine analysis method
if [[ "$1" == "exactTest"* ]]; then
	#TO-DO merge exact test results for selected genotypes
	#Rscript geneSets_vennDiagram.r "$1"
	Rscript exactTest_vennDiagram.r "$inputPath"/"$1"
elif [[ "$1" == "glm"* ]]; then
	Rscript glm_vennDiagram.r "$inputPath"/"$1"
fi

#Rename and move produced plots
for p in *.jpg; do mv $p "$inputPath"; done
