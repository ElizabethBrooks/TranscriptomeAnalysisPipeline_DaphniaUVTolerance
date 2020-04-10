#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash generateBarPlots_binnedCompare_driver.sh twoAlingmentSummaryFiles 
#Usage Ex: bash generateBarPlots_binnedCompare_driver.sh alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv

#Load module necessary for crc servers
module load bio/R/
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve input alignment summary absolute path
inOutputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
inOutDir="$inOutputsPath"/AlignmentsAnalyzed
echo "Plotting $1 and $2 alignment summaries..."
#Plot alignment data using binned bar plots
Rscript alignmentSummary_barPlot_binnedCompare.r "$inOutDir"/"$1" "$inOutDir"/"$2"
echo "Alignment summaries for $1 and $2 have been plotted!"
#Rename produced pdf of plots
outFile1=$(echo "$1" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
outFile2=$(echo "$2" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
outFile=$(echo "$outFile1"_"$outFile2" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
mv "$inOutDir"/plotOverallPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_overallPercentages.jpg
mv "$inOutDir"/plotConcordantPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_concordantPercentages.jpg