#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N stats_edgeR_jobOutput
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash statistics_exactTest_subset.sh countsFile startColPos endColPos
#Usage Ex: bash statistics_exactTest_subset.sh GeneCountAnalysis_subset_run1/geneCounts_merged_counted_htseqTophat2_run1_subset_cleaned.csv 1 6

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
outputCounts="$outputsPath"/geneCountStats_cols"$2"to"$3"
#Retrieve analysis inputs path
inputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
outFile=$(basename "$inputsPath"/"$1" | sed 's/\.csv//g')
#Create directory for output files
mkdir "$outputCounts"_"$outFile"

#Perform DE analysis using edgeR and output analysis results to a txt file
Rscript statistics_exactTest_edgeR.r "$inputsPath"/"$1" $2 $3 > "$outputCounts"_"$outFile"/analysisResults.txt
#Rename and move produced filtered table
mv stats_exactTest.csv "$outputCounts"_"$outFile"/stats_exactTest.csv
#Rename and move produced plots
mv plotBarsBefore.jpg "$outputCounts"_"$outFile"/plotBarsBefore.jpg
mv plotMDSBefore.jpg "$outputCounts"_"$outFile"/plotMDSBefore.jpg
mv plotHeatMapBefore.jpg "$outputCounts"_"$outFile"/plotHeatMapBefore.jpg
mv plotBarsAfter.jpg "$outputCounts"_"$outFile"/plotBarsAfter.jpg
mv plotMDSAfter.jpg "$outputCounts"_"$outFile"/plotMDSAfter.jpg
mv plotHeatMapAfter.jpg "$outputCounts"_"$outFile"/plotHeatMapAfter.jpg
mv plotBCV.jpg "$outputCounts"_"$outFile"/plotBCV.jpg
mv plotMD.jpg "$outputCounts"_"$outFile"/plotMD.jpg
mv plotVolcano.jpg "$outputCounts"_"$outFile"/plotVolcano.jpg
mv plotMA.jpg "$outputCounts"_"$outFile"/plotMA.jpg
