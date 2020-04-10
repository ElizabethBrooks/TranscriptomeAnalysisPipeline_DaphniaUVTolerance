#!/bin/bash
#Usage: Rscript generateBarPlots_binnedMerged_driver.r percentsFile.csv
#Usage Ex: Rscript generateBarPlots_binnedMerged_driver.r alignmentSummarized_legacyTophat2Hisat2_differences_merged.csv

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Move to outputs directory
cd "$outputsPath"/AlignmentsAnalyzed
#Create directory for alignment analysis
outputAnalysis=AlignmentAnalysis
mkdir "$outputAnalysis"
#Script to run Rscripts that generate kMeans plots
if [[ "$1" == *"differences" ]]; then
	#Plot merged data kMeans clustering
	Rscript alignmentSummary_differences_barPlot_binned.r "$1"
else
	#Plot merged data kMeans clustering
	Rscript alignmentSummary_barPlot_binned.r "$1"
fi
#Rename produced plot
outFile=$(echo "$1" | sed 's/\.csv//')
mv Rplots.pdf "$outputAnalysis"/"$outFile".pdf