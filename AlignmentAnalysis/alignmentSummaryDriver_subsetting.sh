#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummaryDriver_subsetting.sh alignmentFolder optionalAssemblyFolder genotypeList optionalRunList
#Usage Ex: bash alignmentSummaryDriver_subsetting.sh aligned_hisat2 PA42 run1 run2
#Usage Ex: bash alignmentSummaryDriver_subsetting.sh aligned_tophat2 PA42 run1 run2 run3
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh trimmed_run1 aligned_hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsHisat2_run1 aligned_hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsTophat2_run1 aligned_hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra

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
#Determine which analysis folder was input
if [[ "$1"  == aligned* ]]; then
	#Retrieve input alignment summary absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	fileName="$1"
	analysisInput=""
elif [[ "$1"  == sorted* || "$1"  == trimmed* ]]; then
	#Retrieve reads input absolute path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve directory name from input folder path
	fileName="$2"
	analysisInput="$1"
else
	echo "ERROR: The input folder of aligned or assembled files were not found... exiting"
	exit 1
fi
#Initialize variables
counter=1
#Loop through all input sets of treatments and perform selected analsysis
for i in "$@"; do
	#Skip first two arguments
	if [ $counter -ge 3 ]; then
		#Determine what type of data folder was input
		if [[ "$1" == trimmed* ]]; then
			inputFolder="$inputsPath"/"$1""$i"_assemblyTrinity/"$fileName"
			genotype="$i"
		elif [[ "$1" == sorted* ]]; then
			inputFolder="$inputsPath"/"$1""$i"_assemblyGenomeTrinity/"$fileName"
			genotype="$i"
		elif [[ "$1" == aligned* ]]; then
			inputFolder="$inputsPath"/"$1"_"$i"
			genotype="$2"
			fileName="$1"_"$i"
		else
			echo "ERROR: Input folder for analysis is not a valid option... exiting!"
			exit 1
		fi
		#Determine what analysis method was used for the input folder of data
		if [[ "$fileName" == *"hisat2"*  ]]; then
			#Set analysis method for folder naming
			analysisMethod="_hisat2"
			analysisArg=$analysisInput$genotype$analysisMethod
			#Set output folder name
			outputStats="$outDir"/alignmentSummarized_"$analysisArg"
			#Retrieve run number for input alignment folder
			runNum=$(echo "$fileName" | sed "s/aligned"$analysisMethod"_run//g")
			#Set header of overall summary csv file
			echo "sample,overall,concordant" > "$outputStats"_run"$runNum".csv
		elif [[ "$fileName" == *"tophat2"* ]]; then
			#Set analysis method for folder naming
			analysisMethod="_tophat2"
			analysisArg=$analysisInput$genotype$analysisMethod
			#Set output folder name
			outputStats="$outDir"/alignmentSummarized_"$analysisArg"
			#Retrieve run number for input alignment folder
			runNum=$(echo "$fileName" | sed "s/aligned"$analysisMethod"_run//g")
			#Set header of overall summary csv file
			echo "sample,leftMapped,rightMapped,overall,concordant" > "$outputStats"_run"$runNum".csv
		else
			echo "ERROR: The $fileName folder of summary files is not valid... exiting"
			exit 1
		fi
		echo "Merging $inputFolder alignment summaries..."
		#Retrieve summaries for each aligned sample
		for f1 in "$inputFolder"/*/; do
			#Retrieve sample name
			sampleName=$(basename "$f1")
			#Retrieve sample summary based on alignment method
			bash alignmentSummary"$analysisMethod"_sample.sh "$f1" "$analysisArg" "$runNum"
			#Combine summaries into one csv file
			cat "$outputStats"_combined_run"$runNum".csv >> "$outputStats"_run"$runNum".csv
			rm "$outputStats"_combined_run"$runNum".csv
		done
		echo "Alignment summaries have been merged!"
		echo "Formatting $inputFolder merged alignment summary..."
		#Run alignment summary formatting
		bash alignmentSummary_formatting.sh "$analysisArg" "$runNum" "$genotype"
		bash alignmentSummary_genotypeMedians.sh "$analysisArg" "$runNum" $numGenotypes "$genotype"
		echo "Merged alignment summary has been formatted!"
	fi
	counter=$(($counter+1))
done
