#!/bin/bash

# BASH script to drive R summary plotting scripts for WGCNA network module analysis

# usage: bash moduleAnalysis_effectsSubset_driver.sh analysisType set
# default usage ex: bash moduleAnalysis_effectsSubset_driver.sh Tolerance OLYM
# usage ex: bash moduleAnalysis_effectsSubset_driver.sh Tolerance Tol
# usage ex: bash moduleAnalysis_effectsSubset_driver.sh Tolerance NTol

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# set DE analysis directory
#Create directory for output files
deDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
deDir=$deDir"/"$analysisType

# retrieve set
set=$2

# minimum module size
minModSize=30

# determine the percent of each module that is associated with each ANOVA effect
Rscript stackedBarPlot_effectSubsets.R $inDir $deDir $set $minModSize
