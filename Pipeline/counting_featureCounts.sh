#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N counting_featureCounts_jobOutput
#$ -pe smp 8
#Script to perform freatureCounts counting of trimmed, aligned, then name sorted
# paired end reads
#Usage: qsub counting_featureCounts.sh sortedNameFolder
#Usage Ex: qsub counting_featureCounts.sh sortedName_samtoolsTophat2_run1
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
analysisTag=".sorted.bam"
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  == sortedName* ]]; then
	#Set name sorted flag (default) with file type flag
	#flags="-f bam"
	sortType="Name"
elif [[ "$1"  == sortedCoordinate* ]]; then
	#Set coordinate sorted flag with file type flag
	#flags="-f bam -r pos"
	sortType="coordinate"
else
	echo "ERROR: The "$1" folder of name or coordinate sorted files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder or bam files were not found... exiting"
	exit 1
fi
#Retrieve sorted reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "counting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/counting://g")
#Move to outputs directory
cd "$outputsPath"
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	outputFolder=counted"$sortType"_featureCounts"$analysisMethod"_run"$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of featureCounts counting of "$1" data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Loop through all sorted forward and reverse paired reads and store the file locations in an array
for f2 in "$inputsPath"/"$1"/*/; do
	#Name of aligned file
	curAlignedSample="$f2"accepted_hits.bam
	#Trim file path from current file name
	curSampleNoPath=$(basename $f2)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/\.bam//')
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Count reads using featureCounts
	echo "Sample $curSampleNoPath is being counted..."
	featureCounts -T 8 -p -t gene -g ID -a "$genomeFile" -s 2 --donotsort -o "$outputFolder"/"$curSampleNoPath"/counted.sam -M "$curAlignedSample"
	echo "Sample $curSampleNoPath has been counted!"
	#Add run inputs to output summary file
	echo "$curSampleNoPath" >> $inputOutFile
	echo featureCounts -T 8 -p -t gene -g ID -a "$genomeFile" -s 2 --donotsort -o "$outputFolder"/"$curSampleNoPath"/counted.sam -M "$curAlignedSample" >> $inputOutFile
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"
