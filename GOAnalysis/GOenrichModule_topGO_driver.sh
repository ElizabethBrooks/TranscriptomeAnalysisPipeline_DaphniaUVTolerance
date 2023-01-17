#!/bin/bash

# BASH script to drive R summary plotting scripts for WGCNA network module analysis

# usage: bash GOenrichModule_topGO_driver.sh analysisType set
# default usage ex: bash GOenrichModule_topGO_driver.sh Tolerance OLYM
# usage ex: bash GOenrichModule_topGO_driver.sh Tolerance Tol
# usage ex: bash GOenrichModule_topGO_driver.sh Tolerance NTol

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# set DE analysis directory
deDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
deDir=$deDir"/"$analysisType

# retrieve functional annotations
GOmaps=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")
GOmaps=$(dirname $GOmaps)
GOmaps=$GOmaps"/geneToGO_tagged_map.txt"

# retrieve set
set=$2

# minimum module size
minModSize=30

# name outputs directory
outDir=$inDir"/GOAnalysis"

# create outputs directory
mkdir $outDir

# determine the direction of expression for each module and effect set
# interaction
Rscript GOenrichModule_topGO_interSubset.R $outDir $deDir $set $minModSize "interaction" $inDir $GOmaps
# treatment
Rscript GOenrichModule_topGO_interSubset.R $outDir $deDir $set $minModSize "treatment" $inDir $GOmaps
# tolerance
Rscript GOenrichModule_topGO_interSubset.R $outDir $deDir $set $minModSize "tolerance" $inDir $GOmaps
