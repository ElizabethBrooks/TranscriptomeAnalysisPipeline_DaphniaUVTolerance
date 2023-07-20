#!/bin/bash

# BASH script to drive GO analysis for DE results

# usage: bash dotplotDE_driver.sh analysisType
# usage ex: bash dotplotDE_driver.sh Tolerance
# default usage ex: bash dotplotDE_driver.sh Genotypes

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")
inDir=$inDir"/"$analysisType

# set outputs directory
outDir=$inDir"/GOAnalysis"

# create outputs directory
mkdir $outDir

# set input paths
positiveTable="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests/"$analysisType"/fisherTest_positiveSelection_modules.csv"

# create a dot plot of significant GO terms for each effect set
Rscript dotplotDE_onlySig.R $outDir $positiveTable
