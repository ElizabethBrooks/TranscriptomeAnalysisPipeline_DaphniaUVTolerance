#!/bin/bash

# script to run Rscripts that log2 transformation of normalized gene counts
# usage: bash logTransformNormalizedCounts_driver.sh analysisType
# usage Ex: bash logTransformNormalizedCounts_driver.sh Genotypes
# usage Ex: bash logTransformNormalizedCounts_driver.sh Tolerance

#Load module for R
#module load bio

# retrieve analysis type
analysisType=$1

#Create directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

# name output directory for the analysis type
inDir=$inDir"/"$analysisType

# run R script to log transform nomalized counts
Rscript logTransformNormalizedCounts.r $inDir
