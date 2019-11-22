#!/bin/bash
#Usage: bash generatePlots_kMeans.sh mergedCounts_file_transposed.csv k
#Usage Ex: bash generatePlots_kMeans.sh mergedCounts_legacy_transposed.csv 3
#Script to run Rscripts that generate kMeans plots
#Retrieve outputs absolute path
outputsFile="TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/outputsPath.txt"
outputsPath=$(head -n 1 $outputsFile)
#Move to outputs directory
cd "$outputsPath"
#Create directory for gene count analysis
outputAnalysis=GeneCountAnalysis
mkdir "$outputAnalysis"
#Plot merged data kMeans clustering
Rscript geneCounts_kMeans.r $1 $2
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf "$outputAnalysis"/"$outFile"_k"$2".pdf