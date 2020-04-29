#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash generateBarPlotsDriver_binnedCompare.sh twoAlingmentSummaryFiles 
#Usage Ex: bash generateBarPlotsDriver_binnedCompare.sh alignmentSummarized_hisat2_run1_formatted.csv alignmentSummarized_tophat2_run2_formatted.csv

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
#Determine what software was used for inputs
if [[ "$1" == *"tophat2"* && "$2" == *"tophat2"* ]]; then
	#Plot alignment data using min intron length binned bar plots
	Rscript alignmentSummary_barPlot_binnedCompareRun.r "$inOutDir"/"$1" "$inOutDir"/"$2"
elif [[ "$1" == *"hisat2"* && "$2" == *"hisat2"* ]]; then
	#Plot alignment data using min intron length binned bar plots
	Rscript alignmentSummary_barPlot_binnedCompareRun.r "$inOutDir"/"$1" "$inOutDir"/"$2"
else
	#Plot alignment data using software binned bar plots
	Rscript alignmentSummary_barPlot_binnedCompareSoftware.r "$inOutDir"/"$1" "$inOutDir"/"$2"
fi
echo "Alignment summaries for $1 and $2 have been plotted!"
#Rename produced pdf of plots
outFile1=$(echo "$1" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
outFile2=$(echo "$2" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
outFile=$(echo "$outFile1"_"$outFile2" | sed 's/_formatted\.csv//' | sed 's/alignmentSummarized_//')
mv "$inOutDir"/plotOverallPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_overallPercentages.jpg
mv "$inOutDir"/plotConcordantPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_concordantPercentages.jpg
