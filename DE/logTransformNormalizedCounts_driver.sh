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

# clean produced tables
file=$inDir"/glmQLF_normalizedCounts_logTransformed.csv"

# fix header
tail -n+2 "$file" > $inDir"/tmpTail.csv"
head -1 "$file" | sed -e 's/^/gene,/' > $inDir"/tmpHeader.csv"

# update table
cat $inDir"/tmpHeader.csv" > "$file"
cat $inDir"/tmpTail.csv" >> "$file"

# clean up
rm $inDir"/tmpHeader.csv"
rm $inDir"/tmpTail.csv"
