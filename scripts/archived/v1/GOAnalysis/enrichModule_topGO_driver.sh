#!/bin/bash

# BASH script to drive GO analysis for WGCNA network modules

# usage: bash enrichModule_topGO_driver.sh analysisType set
# default usage ex: bash enrichModule_topGO_driver.sh Genotypes OLYM
# usage ex: bash enrichModule_topGO_driver.sh Tolerance OLYM
# usage ex: bash enrichModule_topGO_driver.sh Tolerance Tol
# usage ex: bash enrichModule_topGO_driver.sh Tolerance NTol

# retrieve analysis type
analysisType=$1

# set directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# retrieve functional annotations
GOmaps=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")
GOmaps=$(dirname $GOmaps)
GOmaps=$GOmaps"/geneToGO_tagged_map.txt"

# retrieve table of dNdS results
dNdSTable="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/SelectionTests/Pulex_Olympics_kaksResults.csv"

# retrieve set
set=$2

# minimum module size
minModSize=30

# retrieve test type
testType="ks"

# name outputs directory
outDir=$inDir"/GOAnalysis_"$testType"_"$set"_"$minModSize

# create outputs directory
mkdir $outDir

# determine the direction of expression for each module and effect set
Rscript enrichModule_topGO.R $outDir $set $minModSize $inDir $GOmaps $dNdSTable $testType
