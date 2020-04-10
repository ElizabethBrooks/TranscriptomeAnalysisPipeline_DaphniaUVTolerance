#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash generateBarPlots_binnedCompare_driver.sh twoAlingmentSummaryFiles 
#Usage Ex: bash generateBarPlots_binnedCompare_driver.sh alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve input alignment summary absolute path
inputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
outDir="$inputsPath"/AlignmentsAnalyzed
echo "Plotting $1 and $2 alignment summaries..."
#Plot alignment data using binned bar plots
Rscript alignmentSummary_barPlot_binnedCompare.r "$inputsPath"/"$1" "$inputsPath"/"$2"
echo "Alignment summaries for $1 and $2 have been plotted!"
#Rename produced pdf of plots
outFile=$(echo "$1"_"$2" | sed 's/\_formatted\.csv//' | sed 's/alignmentSummarized\_//')
mv Rplots.pdf "$outputAnalysis"/alignmentSummarized_"$outFile"_barPlots.pdf