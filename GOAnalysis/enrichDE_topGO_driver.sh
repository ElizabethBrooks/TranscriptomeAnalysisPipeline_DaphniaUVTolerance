#!/bin/bash

# BASH script to drive GO analysis for DE results

# usage: bash enrichDE_topGO_driver.sh analysisType
# usage ex: bash enrichDE_topGO_driver.sh Tolerance
# default usage ex: bash enrichDE_topGO_driver.sh Genotypes

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inDir=$inDir"/"$analysisType

# retrieve functional annotations
GOmaps=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")
GOmaps=$(dirname $GOmaps)
GOmaps=$GOmaps"/geneToGO_tagged_map.txt"

# retrieve test type
testType="ks"

# set outputs directory
outDir=$inDir"/GOAnalysis_"$testType

# create outputs directory
mkdir $outDir

# determine the direction of expression for each module and effect set
# interaction
Rscript enrichDE_topGO.R $outDir "interaction" $inDir $GOmaps $testType
# treatment
Rscript enrichDE_topGO.R $outDir "treatment" $inDir $GOmaps $testType
# tolerance
Rscript enrichDE_topGO.R $outDir "tolerance" $inDir $GOmaps $testType
