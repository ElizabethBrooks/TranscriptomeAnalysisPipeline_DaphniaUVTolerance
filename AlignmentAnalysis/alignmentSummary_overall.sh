#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_overall.sh alignmentFolders
#Usage Ex: bash alignmentSummary_overall.sh aligned_hisat2_run1 aligned_hisat2_run2 aligned_tophat2_run1 aligned_tophat2_run2 aligned_tophat2_run3

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == *assembly* ]]; then
	analysisInput="assembly"
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve build transcriptome files absolute path
	buildInputsPath="$assemblyPath"/"$1"
	#Retrieve transcriptome reference absolute path for alignment
	buildFile="$assemblyPath"/"$1"/"Trinity.fasta"
	#Retrieve alignment outputs absolute path
	outputsPath="$assemblyPath"/"$1"
	#Retrieve trimmed run folder name used for assembly
	assemblyFolder=$(echo $1 | sed 's/trimmed_run.//')
	trimmedFolder=$(echo $1 | sed "s/$assemblyFolder//")
elif [[ "$1"  == trimmed* ]]; then
	analysisInput="trimmed"
	#Retrieve build genome files absolute path
	buildInputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/building://g")
	#Retrieve genome reference absolute path for alignment
	buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Retrieve alignment outputs absolute path
	outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
	trimmedFolder="$1"
else
	echo "ERROR: The input folder of trimmed or assembled files were not found... exiting"
	exit 1
fi
#Retrieve input alignment summary absolute path
#inputsPath=$(grep "aligning:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory
outDir="$outputsPath"/AlignmentsAnalyzed
#Retrieve folders to analyze from the input arguments to the script
for f1 in $@; do
	#Determine if the folder name was input in the correct format
	if [[ "$f1" == *\/* ]] || [[ "$f1" == *\\* ]]; then
		echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
		exit 1
	fi
	#Determine what analysis method was used for the input folder of data
	if [[ "$f1" == *"hisat2"*  ]]; then
		#Set analysis method for folder naming
		analysisMethod="hisat2"
		#Set output folder name
		outputStats="$outDir"/alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_run//g")
		#Set header of overall summary csv file
		echo "sample,overall,concordant" > "$outputStats"_run"$runNum".csv
	elif [[ "$f1" == *"tophat2"* ]]; then
		#Set analysis method for folder naming
		analysisMethod="tophat2"
		#Set output folder name
		outputStats="$outDir"/alignmentSummarized_"$analysisMethod"
		#Retrieve run number for input alignment folder
		runNum=$(echo "$f1" | sed "s/aligned_"$analysisMethod"_run//g")
		#Set header of overall summary csv file
		echo "sample,leftMapped,rightMapped,overall,concordant" > "$outputStats"_run"$runNum".csv
	else
		echo "ERROR: The $f1 folder or bam files were not found... exiting"
		exit 1
	fi
	echo "Merging $f1 alignment summaries..."
	#Retrieve summaries for each aligned sample
	for f2 in "$inputsPath"/"$f1"/*/; do
		#Retrieve sample name
		sampleName=$(basename "$f2")
		#Retrieve sample summary based on alignment method
		bash alignmentSummary_"$analysisMethod"_sample.sh "$f2" "$analysisMethod" "$runNum"
		#Combine summaries into one csv file
		cat "$outputStats"_combined_run"$runNum".csv >> "$outputStats"_run"$runNum".csv
		rm "$outputStats"_combined_run"$runNum".csv
	done
	echo "Alignment summaries for $f1 have been merged!"
	echo "Formatting $f1 merged alignment summary..."
	#Run alignment summary formatting
	bash alignmentSummary_formatting.sh "$analysisMethod" "$runNum"
	echo "Merged alignment summary for $f1 has been formatted!"
done