#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N sorting_samtools_jobOutput
#$ -pe smp 8
#Script to perform samtools sorting of trimmed, then aligned
# paired end reads
#Usage: qsub sorting_samtools.sh sortingTarget sortingMethod alignedFolder optionalAssembledFolder
#Usage Ex: qsub sorting_samtools.sh genome name aligned_tophat2_run1
#Usage Ex: qsub sorting_samtools.sh assembly name aligned_hisat2_run1 sortedCoordinate_samtoolsTophat2_run1E05_assemblyGenomeTrinity

#Required modules for ND CRC servers
module load bio/2.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve sorting method flags from input
if [[ "$2" == "name" || "$2" == "Name" || "$2" == "n" || "$2" == "N" ]]; then
	#Name sorted flag with num threads flag
	flags="-@ 8 -n"
	methodTag="Name"
elif [[ "$2" == "coordinate" || "$2" == "Coordinate" || "$2" == "c" || "$2" == "C" ]]; then
	#Coordinate sorted with num threads flag
	flags="-@ 8"
	methodTag="Coordinate"
else
	#Report error with input flag
	echo "ERROR: a flag for sorting method (name or coordiante) is expected... exiting"
	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == assembly ]]; then
	#inputsPath reads input absolute path
	inputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	inputsPath="$inputsPath"/"$4"
	#Retrieve alignment outputs absolute path
	outputsPath="$inputsPath"
elif [[ "$1"  == genome ]]; then
	#Retrieve aligned reads input absolute path
	inputsPath=$(grep "aligning:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
	#Retrieve sorting outputs absolute path
	outputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
else
	echo "ERROR: The input alignment target is not valid... exiting!"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$3" == *"hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$3" == *"tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The "$3" folder of files were not found... exiting"
	exit 1
fi
#Move to outputs directory
cd "$outputsPath"
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
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
		echo "Creating folder for run $runNum of Samtools sorting of "$3" data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Loop through all reads and sort sam/bam files for input to samtools
for f1 in "$inputsPath"/"$3"/*/; do
	#Determine what extension the files have
	#curSampleHits=$(echo "$f1"*)
	#curSampleHits=$(basename "$curSampleHits")
	#extension=${curSampleHits##*.}
	#Name of aligned file
	#curAlignedSample="$f1"accepted_hits."$extension"
	curAlignedSample="$f1"accepted_hits.bam
	#Trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Output current sample name to summary file
	echo "$curSampleNoPath" >> $inputOutFile
	#Run samtools to prepare mapped reads for sorting by name
	#using 8 threads
	echo "Sample $curSampleNoPath is being name sorted..."
	samtools sort -@ 8 -n -o "$outputFolder"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample"
	echo "Sample $curSampleNoPath has been name sorted!"
	#Add run inputs to output summary file
	echo samtools sort -@ 8 -n -o "$outputFolder"/"$curSampleNoPath"/sortedName.bam -T /tmp/"$curSampleNoPath".sortedName.bam "$curAlignedSample" >> "$inputOutFile"
	#Determine which sorting method is to be performed
	if [[ "$methodTag" == "Coordinate" ]]; then
		#Run fixmate to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		#Clean up
		rm "$outputFolder"/"$curSampleNoPath"/sortedName.bam
		#Run samtools to prepare mapped reads for sorting by coordinate
		#using 8 threads
		echo "Sample $curSampleNoPath is being sorted..."
		samtools sort "$flags" -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam.bam -T /tmp/"$curSampleNoPath".sorted.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		echo "Sample $curSampleNoPath has been sorted!"
		rm "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam
		#Add run inputs to output summary file
		echo samtools fixmate "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam >> "$inputOutFile"
		echo samtools sort "$flags" -o "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam.bam -T /tmp/"$curSampleNoPath".sorted.bam "$outputFolder"/"$curSampleNoPath"/sortedFixed.bam >> "$inputOutFile"
	else
		#Run fixmate to update paired-end flags for singletons
		echo "Sample $curSampleNoPath singleton flags are being updated..."
		samtools fixmate "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam.bam
		echo "Sample $curSampleNoPath singleton flags have been updated!"
		rm "$outputFolder"/"$curSampleNoPath"/sortedName.bam
		#Add run inputs to output summary file
		echo samtools fixmate "$outputFolder"/"$curSampleNoPath"/sortedName.bam "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam.bam >> "$inputOutFile"
	fi
done
#Copy previous summaries
cp "$inputsPath"/"$3"/*.txt "$outputFolder"
