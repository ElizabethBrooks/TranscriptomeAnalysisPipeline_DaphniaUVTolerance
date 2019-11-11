#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N build_hisat2_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/hisat2/2.1.0
#Prepare for alignment
cd ..
dirFlag=0
runNum=1
buildFile=$(tail -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt")
#Build reference genome if folder does not exist
#Build output directory
outputFolder="reference_hisat2_build"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Build files already exsist
	echo "Build files for hisat2 already found in $outputFolder... exiting"
	exit 1
else
	#Build files do not exsist
	echo "Creating $outputFolder folder for hisat2 build..."
fi
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
if [ $? -eq 0 ]; then
	#Trim file path from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNoPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g')
	#Copy genome build fasta file to hisat2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNoPath"
	#Begin hisat2 build
	echo "Beginning hisat2 build... "
	hisat2-build -p 8 -f "$outputFolder"/"$buildFileNoPath" "$outputFolder"/"$buildFileNoPath"
	echo "hisat2 build complete!"
else
	echo "Build folder reference_hisat2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo hisat2-build -p 8 -f "$outputFolder"/"$buildFileNoPath" "$outputFolder"/"$buildFileNoPath" >> $inputOutFile