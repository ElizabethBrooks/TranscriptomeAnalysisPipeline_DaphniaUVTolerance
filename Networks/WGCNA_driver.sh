#!/bin/bash
#Script to run Rscripts that perform WGCNA analysis of OLYM data
#Usage: bash WGCNA_driver_$analysisType.sh
#Usage ex: bash WGCNA_driver.sh Tolerance
#Usage ex: bash WGCNA_driver.sh Genotypes

# retrieve analysis type
analysisType=$1

# set the soft-thresholding power
# it is recommended for less than 20 samples to use
# 9 for unsigned and 18 for signed
softThresh=18

# set the minimum module size
minModSize=30

#Create directory for output files
inDir=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

# name output directory for the analysis type
inDir=$inDir"/"$analysisType

# input counts file
inputCounts=$inDir"/glmQLF_normalizedCounts_logTransformed.csv"

# get current directory
currDir=$(pwd)

# retrieve experimental design data path
cd $(dirname "../InputData/expDesign_WGCNA_Olympics.csv")
#cd $(dirname "../InputData/expDesign_treatment_WGCNA_Olympics.csv")
designPath=$(pwd)
expDesign=$designPath"/expDesign_WGCNA_Olympics.csv"
#expDesign=$designPath"/expDesign_treatment_WGCNA_Olympics.csv"

# move back to scripts directory
cd $currDir

# name output directories
outPath=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
outDir=$outPath"/"$analysisType

# create output directories
mkdir $outPath
mkdir $outDir

# WGCNA data input
Rscript "dataInput_consensus"$analysisType"_WGCNA.R" $outDir $inputCounts 1 24 $expDesign $analysisType
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 OLYM $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 12 Tol $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 13 24 NTol $expDesign
#Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 6 Y05 $expDesign
#Rscript dataInput_set_WGCNA.R $outDir $inputCounts 7 12 Y023 $expDesign
#Rscript dataInput_set_WGCNA.R $outDir $inputCounts 13 18 E05 $expDesign
#Rscript dataInput_set_WGCNA.R $outDir $inputCounts 19 24 R2 $expDesign

# WGCNA soft power picking
Rscript pickSoftPower_consensus_WGCNA.R $outDir $analysisType
Rscript pickSoftPower_set_WGCNA.R $outDir OLYM
Rscript pickSoftPower_set_WGCNA.R $outDir Tol
Rscript pickSoftPower_set_WGCNA.R $outDir NTol
#Rscript pickSoftPower_set_WGCNA.R $outDir Y05
#Rscript pickSoftPower_set_WGCNA.R $outDir Y023
#Rscript pickSoftPower_set_WGCNA.R $outDir E05
#Rscript pickSoftPower_set_WGCNA.R $outDir R2

# WGCNA network construction
Rscript "networkConstruction_consensus"$analysisType"_WGCNA.R" $outDir $softThresh $minModSize
Rscript networkConstruction_set_WGCNA.R $outDir OLYM 8 $minModSize
Rscript networkConstruction_set_WGCNA.R $outDir Tol 14 $minModSize
Rscript networkConstruction_set_WGCNA.R $outDir NTol 14 $minModSize
#Rscript networkConstruction_set_WGCNA.R $outDir Y05 20 $minModSize
#Rscript networkConstruction_set_WGCNA.R $outDir Y023 12 $minModSize
#Rscript networkConstruction_set_WGCNA.R $outDir E05 14 $minModSize
#Rscript networkConstruction_set_WGCNA.R $outDir R2 9 $minModSize

# WGCNA consensus network analsysis
Rscript networkAnalysis_consensus_WGCNA.R $outDir OLYM $minModSize $analysisType
Rscript networkAnalysis_consensus_WGCNA.R $outDir Tol $minModSize $analysisType
Rscript networkAnalysis_consensus_WGCNA.R $outDir NTol $minModSize $analysisType
#Rscript networkAnalysis_consensus_WGCNA.R $outDir Y05 $minModSize $analysisType
#Rscript networkAnalysis_consensus_WGCNA.R $outDir Y023 $minModSize $analysisType
#Rscript networkAnalysis_consensus_WGCNA.R $outDir E05 $minModSize $analysisType
#Rscript networkAnalysis_consensus_WGCNA.R $outDir R2 $minModSize $analysisType

# WGCNA compare consensus eigengene networks
Rscript "eigengeneNetworks_consensus"$analysisType"_WGCNA.R" $outDir $minModSize

# WGCNA export eigengene expression values
Rscript eigengeneExpression_set_WGCNA.R $outDir OLYM $minModSize
Rscript eigengeneExpression_set_WGCNA.R $outDir Tol $minModSize
Rscript eigengeneExpression_set_WGCNA.R $outDir NTol $minModSize
#Rscript eigengeneExpression_set_WGCNA.R $outDir Y05 $minModSize
#Rscript eigengeneExpression_set_WGCNA.R $outDir Y023 $minModSize
#Rscript eigengeneExpression_set_WGCNA.R $outDir E05 $minModSize
#Rscript eigengeneExpression_set_WGCNA.R $outDir R2 $minModSize
