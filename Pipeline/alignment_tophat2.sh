#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2_jobOutput
#$ -pe smp 8
#Script to perform tophat2 alignment of trimmed
# paired end reads
#Note that a bowtie2 genome refernce build folder needs to be generated first
#Usage: qsub alignment_tophat2.sh trimmedFolder minIntronLength maxIntronLength
#Usage Ex: qsub alignment_tophat2.sh trimmed_run1 4 14239
#Default usage Ex: qsub alignment_tophat2.sh trimmed_run1

#Required modules for ND CRC servers
module load bio
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
if [[ "$1"  != trimmed* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Retrieve trimmed reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
#Retrieve build genome files absolute path
buildInputsPath=$(grep "building:" ../InputData/outputPaths.txt | tr -d " " | sed "s/building://g")
#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "aligning:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
#Move to outputs directory
cd "$outputsPath"
#Prepare for alignment
dirFlag=0
runNum=1
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	#Tophat output directory name
	outputFolder="aligned_tophat2_run$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of tophat2 alignment on "$1" data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Build output directory for Tophat reference
buildOut="$buildInputsPath"/"reference_bowtie2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/')
buildFileNoEx=$(echo $buildFileNoEx | sed 's/\.fa//')
#Loop through all forward and reverse paired reads and run tophat2 on each pair
# using 8 threads
for f1 in "$inputsPath"/"$1"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoEx=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Determine if intron lengths were entered
	if [[ -z "$2" || -z "$3" ]]; then #Arguments were not entered
		#Run tophat2 with default settings
		echo "Sample $curSampleNoEx is being aligned..."
		tophat2 -p 8 -G "$genomeFile" -o "$outputFolder"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "tophat2 -p 8 -G" "$genomeFile" -o "$outputFolder"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz >> $inputOutFile
	else #Run tophat2 using input intron lengths
		echo "Sample $curSampleNoEx is being aligned..."
		tophat2 -p 8 -i $2 -I $3 -G "$genomeFile" -o "$outputFolder"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "tophat2 -p 8 -i $2 -I $3 -G" "$genomeFile" -o "$outputFolder"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz >> $inputOutFile
	fi
	echo "Sample $curSampleNoEx has been aligned!"
done
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"
cp "$buildOut"/*.txt "$outputFolder"
