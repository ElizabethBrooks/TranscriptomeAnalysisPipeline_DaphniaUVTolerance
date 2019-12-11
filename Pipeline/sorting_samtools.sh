#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sorting_samtools_jobOutput
#$ -pe smp 8
#Script to perform samtools sorting of trimmed, then aligned
# paired end reads
#Usage: qsub sorting_samtools.sh -sortingMethod alignedFolder
#Usage Ex: qsub sorting_samtools.sh -name aligned_tophat2_run1
#Required modules for ND CRC servers
module load bio
#module load bio/python/2.7.14
#module load bio/htseq/0.11.2
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve sorting method flags from input
if [[ "$1" == "-name" || "$1" == "-Name" || "$1" == "-n" || "$1" == "-N" ]]; then
	#Name sorted flag with num threads flag
	flags="-@ 8 -n"
	methodTag="Name"
elif [[ "$1" == "-coordinate" || "$1" == "-Coordinate" || "$1" == "-c" || "$1" == "-C" ]]; then
	#Coordinate sorted with num threads flag
	flags="-@ 8"
	methodTag="Coordinate"
else
	#Report error with input flag
	echo "ERROR: a flag for sorting method (name or coordiante) is expected... exiting"
	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$2" == *\/* ]] || [[ "$2" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$2"  != aligned* ]]; then
	echo "ERROR: The "$2" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$2" == *"hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$2" == *"tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The "$2" folder of files were not found... exiting"
	exit 1
fi
#Retrieve aligned reads input absolute path
inputsPath=$(grep "alignment:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignment://g")
#Retrieve sorting outputs absolute path
outputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Move to outputs directory
cd "$outputsPath"
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	outputFolder=sorted"$methodTag"_samtools"$analysisMethod"_run"$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of Samtools sorting of "$2" data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Sort input sam/bam files if folder does not already exist
if [ $? -eq 0 ]; then
	echo "Creating folder for sorted files..."
	#Loop through all reads and sort sam/bam files for input to samtools
	for f1 in "$inputsPath"/"$2"/*/; do
		#Determine what extension the files have
		curSampleHits=$(echo "$f1"*)
		curSampleHits=$(basename "$curSampleHits")
		extension=${curSampleHits##*.}
		#Name of aligned file
		curAlignedSample="$f1"accepted_hits."$extension"
		#Trim file path from current folder name
		curSampleNoPath=$(echo "$f1")
		#Create directory for current sample outputs
		mkdir "$outputFolder"/"$curSampleNoPath"
		#Run samtools to prepare mapped reads for sorting by name
		#using 8 threads
		echo "Sample $curSampleNoPath is being sorted..."
		samtools sort "$flags" -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$curAlignedSample"
		echo "Sample $curSampleNoPath has been sorted!"
		#Add run inputs to output summary file
		echo "$curSampleNoPath" >> $inputOutFile
		echo samtools sort "$flags" -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam -T /tmp/"$curSampleNoPath".sorted.bam "$curAlignedSample" >> $inputOutFile
	done
	#Copy previous summaries
	cp "$inputsPath"/"$2"/*.txt "$outputFolder"
else
	echo "Sorted files already exists, skipping sorting..."
fi
