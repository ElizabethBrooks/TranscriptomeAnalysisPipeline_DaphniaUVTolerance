#!/bin/bash

# BASH script to drive GO analysis for DE results

# usage: bash dotplot_driver.sh analysisType
# default usage ex: bash dotplot_driver.sh Tolerance

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inDir=$inDir"/"$analysisType

# set outputs directory
outDir=$inDir"/GOAnalysis"

# create outputs directory
mkdir $outDir

# create a dot plot of significant GO terms for each effect set
Rscript dotplot_onlySig.R $outDir