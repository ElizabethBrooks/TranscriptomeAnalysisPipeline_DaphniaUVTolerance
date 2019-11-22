#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N build_bowtie2_jobOutput
#Required modules for ND CRC servers
module load bio
#Prepare for alignment
dirFlag=0
runNum=1
#Retrieve genome file path for building
buildFile=$(tail -n 1 "InputData/genomeFilePaths.txt")
#Retrieve outputs absolute path
outputsFile="InputData/outputsPath.txt"
outputsPath=$(head -n 1 $outputsFile)
#Move to outputs directory
cd "$outputsPath"
#Create output directory
outputFolder="reference_bowtie2_build"
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Build files already exsist
	echo "Build files for bowtie2 already found in $outputFolder... exiting"
	exit 1
else
	#Build files do not exsist
	echo "Creating $outputFolder folder for bowtie2 build..."
fi
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
if [ $? -eq 0 ]; then
	#Trim file path from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNoPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g')
	#Copy genome build fasta file to bowtie2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNoPath"
	#Begin Bowtie2 build
	echo "Beginning bowtie2 build... "
	bowtie2-build "$outputFolder"/"$buildFileNoPath" "$outputFolder"/"$buildFileNoPath"
	echo "Bowtie2 build complete!"
else
	echo "Build folder reference_bowtie2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo bowtie2-build "$outputFolder"/"$buildFileNoPath" "$outputFolder"/"$buildFileNoPath" >> $inputOutFile