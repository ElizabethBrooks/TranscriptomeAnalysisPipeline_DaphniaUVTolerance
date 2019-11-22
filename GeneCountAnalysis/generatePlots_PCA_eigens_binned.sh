#!/bin/bash
#Usage: bash generatePlots_PCA_eigens_binned.sh mergedCounts_file_transposed.csv bin
#Usage Ex: bash generatePlots_PCA_eigens_binned.sh mergedCounts_legacy_transposed.csv treatment
#Script to run Rscripts that generate binned PCA plots
#Retrieve outputs absolute path
outputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/outputsPath.txt"
outputsPath=$(head -n 1 $outputsFile)
#Move to outputs directory
cd "$outputsPath"
#Create directory for gene count analysis
outputAnalysis=GeneCountAnalysis
mkdir "$outputAnalysis"
#Plot merged data binned PCA
Rscript geneCounts_PCA_eigens_binned.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf "$outputAnalysis"/"$outFile"_bin"$2".pdf