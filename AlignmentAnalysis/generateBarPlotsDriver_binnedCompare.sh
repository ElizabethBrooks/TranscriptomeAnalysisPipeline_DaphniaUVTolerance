#!/bin/bash
#Script to generate venn diagrams from edgeR stats
#Usage: bash generateBarPlotsDriver_binnedCompare.sh twoAlingmentSummaryFiles 
#Usage Ex: bash generateBarPlotsDriver_binnedCompare.sh hisat2_run1 tophat2_run2
#Alternate usage Ex: bash generateBarPlotsDriver_binnedCompare.sh trimmed_run1E05_hisat2_run1 hisat2_run1
#Alternate usage Ex: bash generateBarPlotsDriver_binnedCompare.sh sortedCoordinate_samtoolsHisat2_run1E05_hisat2_run1 hisat2_run1
#Alternate usage Ex: bash generateBarPlotsDriver_binnedCompare.sh PA42_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1 hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra
#Alternate usage Ex: bash generateBarPlotsDriver_binnedCompare.sh sortedCoordinate_samtoolsHisat2_run1 hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra

#Load module necessary for crc servers
#module load R
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve input alignment summary absolute path
inOutputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory and file
inOutDir="$inOutputsPath"/AlignmentsAnalyzed
echo "Plotting alignment summaries..."
#Check number of input folders
counter=1
fileList=""
if [ $# -ge 3 ]; then
	#Determine comparison to be performed
	if [[ "$1" == PA42* ]]; then
		#Generate file list
		for i in "$@"; do
			#Skip first 3 arguments
			if [ $counter -ge 4 ]; then
				currFile=alignmentSummarized_"$2$i"_"$3"_medians.csv
				fileList="$fileList$inOutDir/$currFile "
				printList="$fileList$currFile "
			fi
			counter=$((counter+1))
		done
		addAnalysis=$inOutDir/alignmentSummarized_"$1"_medians.csv
		fileList="$fileList$addAnalysis"
		outFile="$1"_vs_"$2"TargetGenotypes_"$3"
	else
		#Generate file list
		for i in "$@"; do
			#Skip first 2 arguments
			if [ $counter -ge 3 ]; then
				currFile=alignmentSummarized_"$1$i"_"$2"_medians.csv
				fileList="$fileList$inOutDir/$currFile "
				printList="$fileList$currFile "
			fi
			counter=$((counter+1))
		done
		outFile="$1"TargetGenotypes_"$2"
	fi
	#Plot alignment data using software binned bar plots
	Rscript alignmentSummary_barPlotMedians_binnedCompareGenotype.r $fileList
	#Rename produced pdf of plots
	mv "$inOutDir"/plotMedianOverallPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_medianOverallPercentages.jpg
	mv "$inOutDir"/plotMedianConcordantPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_medianConcordantPercentages.jpg
else
	printList="alignmentSummarized_"$1" alignmentSummarized_"$2
	#Determine what software was used for inputs
	if [[ "$1" == *"Trinity"* && "$2" != *"Trinity"* ]]; then
		#Plot alignment data using target (assembly vs aligned) binned bar plots
		Rscript alignmentSummary_barPlot_binnedCompareTarget.r "$inOutDir"/alignmentSummarized_"$1"_formatted.csv "$inOutDir"/alignmentSummarized_"$2"_formatted.csv
		Rscript alignmentSummary_barPlotMedians_binnedCompareTarget.r "$inOutDir"/alignmentSummarized_"$1"_medians.csv "$inOutDir"/alignmentSummarized_"$2"_medians.csv
	elif [[ "$1" == *"tophat2"* && "$2" == *"tophat2"* ]]; then
		#Plot alignment data using run binned bar plots
		Rscript alignmentSummary_barPlot_binnedCompareRun.r "$inOutDir"/alignmentSummarized_"$1"_formatted.csv "$inOutDir"/alignmentSummarized_"$2"_formatted.csv
		Rscript alignmentSummary_barPlotMedians_binnedCompareRun.r "$inOutDir"/alignmentSummarized_"$1"_medians.csv "$inOutDir"/alignmentSummarized_"$2"_medians.csv
	elif [[ "$1" == *"hisat2"* && "$2" == *"hisat2"* ]]; then
		#Plot alignment data using run binned bar plots
		Rscript alignmentSummary_barPlot_binnedCompareRun.r "$inOutDir"/alignmentSummarized_"$1"_formatted.csv "$inOutDir"/alignmentSummarized_"$2"_formatted.csv
		Rscript alignmentSummary_barPlotMedians_binnedCompareRun.r "$inOutDir"/alignmentSummarized_"$1"_medians.csv "$inOutDir"/alignmentSummarized_"$2"_medians.csv
	else
		#Plot alignment data using software binned bar plots
		Rscript alignmentSummary_barPlot_binnedCompareSoftware.r "$inOutDir"/alignmentSummarized_"$1"_formatted.csv "$inOutDir"/alignmentSummarized_"$2"_formatted.csv
		Rscript alignmentSummary_barPlotMedians_binnedCompareSoftware.r "$inOutDir"/alignmentSummarized_"$1"_medians.csv "$inOutDir"/alignmentSummarized_"$2"_medians.csv
	fi
	echo "Software binned alignment summaries for $printList have been plotted!"
	#Rename produced pdf of plots
	outFile="$1"_"$2"
	mv "$inOutDir"/plotOverallPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_overallPercentages.jpg
	mv "$inOutDir"/plotConcordantPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_concordantPercentages.jpg
	mv "$inOutDir"/plotMedianOverallPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_medianOverallPercentages.jpg
	mv "$inOutDir"/plotMedianConcordantPercentages.jpg "$inOutDir"/alignmentSummarized_"$outFile"_medianConcordantPercentages.jpg
fi