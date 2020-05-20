#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummaryDriver_subsetting.sh alignmentFolders
#Usage Ex: bash alignmentSummaryDriver_subsetting.sh aligned_hisat2_run1 aligned_hisat2_run2 aligned_tophat2_run1 aligned_tophat2_run2 aligned_tophat2_run3
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh trimmed_run1E05_assemblyTrinity/aligned_hisat2_run1 trimmed_run1R2_assemblyTrinity/aligned_hisat2_run1 trimmed_run1Y023_5_assemblyTrinity/aligned_hisat2_run1 trimmed_run1Y05_assemblyTrinity/aligned_hisat2_run1 trimmed_run1PA_assemblyTrinity/aligned_hisat2_run1 trimmed_run1Sierra_assemblyTrinity/aligned_hisat2_run1
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsHisat2_run1E05_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1R2_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1Y023_5_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1Y05_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1PA_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyGenomeTrinity/aligned_hisat2_run1
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsTophat2_run1E05_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1R2_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1Y023_5_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1Y05_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1PA_assemblyGenomeTrinity/aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1Sierra_assemblyGenomeTrinity/aligned_hisat2_run1

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory
outDir="$outputsPath"/AlignmentsAnalyzed
#Set number of genotypes
numGenotypes=6
#Retrieve folders to analyze from the input arguments to the script
for f1 in $@; do
	#Determine which analysis folder was input
	if [[ "$f1"  == aligned* ]]; then
		#Retrieve input alignment summary absolute path
		inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
		analysisInput=""
	elif [[ "$f1"  == *assembly* ]]; then
		#Retrieve reads input absolute path
		basePath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
		inputsPath=$(dirname "$basePath"/"$f1")
		f1=$(basename "$basePath"/$f1)
		#Retrieve directory name from input folder path
		analysisInput=$(basename "$inputsPath")
		analysisInput=$analysisInput"_"
	else
		echo "ERROR: The input folder of aligned or assembled files were not found... exiting"
		exit 1
	fi
	#Determine what analysis method was used for the input folder of data
	if [[ "$f1" == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
		analysisArg=$analysisInput$analysisMethod
		#Set output folder name
		outputStats="$outDir"/alignmentSummarized_"$analysisArg"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_run//g")
		#Set header of overall summary csv file
		echo "sample,overall,concordant" > "$outputStats"_run"$runNum".csv
	elif [[ "$f1" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
		analysisArg=$analysisInput$analysisMethod
		#Set output folder name
		outputStats="$outDir"/alignmentSummarized_"$analysisArg"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_run//g")
		#Set header of overall summary csv file
		echo "sample,leftMapped,rightMapped,overall,concordant" > "$outputStats"_run"$runNum".csv
	else
		echo "ERROR: The $f1 folder of bam files were not found... exiting"
		exit 1
	fi
	echo "Merging $analysisArg"/"$f1 alignment summaries..."
	#Retrieve summaries for each aligned sample
	for f2 in "$inputsPath"/"$f1"/*/; do
		#Retrieve sample name
		sampleName=$(basename "$f2")
		#Retrieve sample summary based on alignment method
		bash alignmentSummary_"$analysisMethod"_sample.sh "$f2" "$analysisArg" "$runNum"
		#Combine summaries into one csv file
		cat "$outputStats"_combined_run"$runNum".csv >> "$outputStats"_run"$runNum".csv
		rm "$outputStats"_combined_run"$runNum".csv
	done
	echo "Alignment summaries for $analysisArg"/"$f1 have been merged!"
	echo "Formatting $analysisArg"/"$f1 merged alignment summary..."
	#Run alignment summary formatting
	bash alignmentSummary_formatting.sh "$analysisArg" "$runNum"
	bash alignmentSummary_genotypeMedians.sh "$analysisArg" "$runNum" $numGenotypes
	echo "Merged alignment summary for $analysisArg"/"$f1 has been formatted!"
done