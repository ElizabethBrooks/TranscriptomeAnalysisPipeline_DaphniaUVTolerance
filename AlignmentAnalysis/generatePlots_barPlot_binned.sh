#!/bin/bash
#Usage: Rscript generatePlots_barPlot_binned.r percentsFile.csv
#Usage Ex: Rscript generatePlots_barPlot_binned.r alignmentSummarized_legacyTophat2Hisat2_differences_merged.csv
#Script to run Rscripts that generate kMeans plots
if [[ "$1" == *"differences" ]]; then
	#Plot merged data kMeans clustering
	Rscript alignmentSummary_differences_barPlot_binned.r $1
else
	#Plot merged data kMeans clustering
	Rscript alignmentSummary_barPlot_binned.r $1
fi
#Rename produced plot
outFile=$(echo $1 | sed 's/\.csv//')
mv Rplots.pdf "$outFile".pdf