#!/bin/bash

# BASH script to plot GO analysis results for WGCNA network modules

# usage: bash dotplotModule_driver.sh analysisType
# usage ex: bash dotplotModule_driver.sh Tolerance OLYM
# usage ex: bash dotplotModule_driver.sh Tolerance Tol
# usage ex: bash dotplotModule_driver.sh Tolerance NTol

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# retrieve set
set=$2

# minimum module size
minModSize=30

# name outputs directory
outDir=$inDir"/GOAnalysis_"$set"_"$minModSize

# create outputs directory
#mkdir $outDir

# create a dot plot of significant GO terms for each effect set
#Rscript dotplotModule_onlySig.R $outDir $set $minModSize $inDir

# create a dot plot of most significant GO terms for each effect set
Rscript dotplotModule_onlyMostSig.R $outDir $set $minModSize $inDir
