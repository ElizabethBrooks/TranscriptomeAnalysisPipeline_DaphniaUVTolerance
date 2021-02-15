#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummaryDriver_subsetting.sh alignmentFolder optionalAssemblyFolder genotypeList optionalRunList
#Usage Ex: bash alignmentSummaryDriver_subsetting.sh aligned_hisat2 run4
#Usage Ex: bash alignmentSummaryDriver_subsetting.sh aligned_tophat2 run1 run2 run3
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh trimmed_run1 aligned_hisat2_run1 E05 Y05 R2 Y023_5 PA Sierra
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsHisat2_run1 aligned_hisat2_run1 PA42_v4.1 E05 Y05 R2 Y023_5 PA Sierra
#Alternate usage Ex: bash alignmentSummaryDriver_subsetting.sh sortedCoordinate_samtoolsHisat2_run2 aligned_hisat2_run1 PA42_v3.0 E05 Y05 R2 Y023_5 PA Sierra

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == aligned* ]]; then
	#Retrieve input alignment summary absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	fileName="$1"
	analysisInput=""
	#Set number of genotypes
	numGenotypes=6
elif [[ "$1"  == sorted* || "$1"  == trimmed* ]]; then
	if [[ "$1" == trimmed* ]]; then
		#Retrieve reads input absolute path
		inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	elif [[ "$1" == sorted* ]]; then
		#Retrieve reads input absolute path
		inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
		#Set genome tag
		genomeTag="$3"
	fi
	#Retrieve directory name from input folder path
	fileName="$2"
	analysisInput="$1"
	#Set number of genotypes
	numGenotypes=1
else
	echo "ERROR: The input folder of aligned or assembled files were not found... exiting"
	exit 1
fi
#Initialize variables
counter=1
#Loop through all input sets of treatments and perform selected analsysis
analysisFlag=-1
for i in "$@"; do
	if [[ "$1" == aligned* ]]; then
		#Skip first argument
		if [ $counter -ge 2 ]; then
			#Set input data folder
			inputFolder="$inputsPath"/"$1"_"$i"
			genotype="$2"
			fileName="$1"_"$i"
			#Set analysis flag
			analysisFlag=1
		fi
	elif [[ "$1" == trimmed*  ]]; then
		#Skip first two arguments
		if [ $counter -ge 3 ]; then
			#Set input data folder
			inputFolder="$inputsPath"/"$1""$i"_assemblyTrinity/"$fileName"
			genotype="$i"
			#Set analysis flag
			analysisFlag=1
		fi
	else
		#Skip first three arguments
		if [ $counter -ge 4 ]; then
			#Set input data folder
			inputFolder="$inputsPath"/"$1""$i"_assembly"$genomeTag"Trinity/"$fileName"
			genotype="$i"
			#Set analysis flag
			analysisFlag=1
		fi
	fi
	#Determine if analysis should begin
	if [[ "$analysisFlag" == 1 ]]; then
		#Set outputs directory
		outputsPath=$(dirname "$inputFolder")
		outDir="$outputsPath"/AlignmentsAnalyzed
		mkdir "$outDir"
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
			bash alignmentSummary"$analysisMethod"_sample.sh "$f1" "$analysisArg" "$runNum" "$outputStats"
			#Combine summaries into one csv file
			cat "$outputStats"_combined_run"$runNum".csv >> "$outputStats"_run"$runNum".csv
			rm "$outputStats"_combined_run"$runNum".csv
		done
		echo "Alignment summaries have been merged!"
		echo "Formatting merged alignment summary..."
		#Run alignment summary formatting
		bash alignmentSummary_formatting.sh "$analysisArg" "$runNum" "$genotype" "$outputStats"
		echo "Merged alignment summary has been formatted!"
		#Generate median values for each genotype
		echo "Generating median values..."
		bash alignmentSummary_genotypeMedians.sh "$analysisArg" "$runNum" "$numGenotypes" "$outputStats"
		echo "Median values have been generated!"
	fi
	counter=$(($counter+1))
done
