#!/bin/bash

# BASH script to drive GO analysis for DE results

# usage: bash enrichDE_topGO_driver.sh analysisType
# default usage ex: bash enrichDE_topGO_driver.sh Tolerance

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inDir=$inDir"/"$analysisType

# retrieve functional annotations
GOmaps=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")
GOmaps=$(dirname $GOmaps)
GOmaps=$GOmaps"/geneToGO_tagged_map.txt"

# set outputs directory
outDir=$inDir"/GOAnalysis"

# create outputs directory
mkdir $outDir

# determine the direction of expression for each module and effect set
# interaction
Rscript enrichDE_topGO.R $outDir "interaction" $inDir $GOmaps
# treatment
Rscript enrichDE_topGO.R $outDir "treatment" $inDir $GOmaps
# tolerance
Rscript enrichDE_topGO.R $outDir "tolerance" $inDir $GOmaps
