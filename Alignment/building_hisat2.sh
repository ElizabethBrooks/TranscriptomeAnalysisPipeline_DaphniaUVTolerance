#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N building_hisat2_jobOutput
#$ -pe smp 4

#Script to generate a hisat2 genome refernce build folder
#Usage: qsub building_hisat2.sh trimmedFolder
#Usage ex: qsub building_hisat2.sh trimmed
#Usage ex: qsub building_hisat2.sh trimmed_run1

#Required modules for ND CRC servers
module load bio/2.0
#module load bio/hisat2/2.1.0

#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve build outputs absolute path
outputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/buildingGenome://g")

#Move to outputs directory
cd "$outputsPath"
#Create output directory
outputFolder="reference_hisat2_build"
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Add software version to output summary file
hisat2-build --version > $inputOutFile

#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
if [ $? -eq 0 ]; then
	#Trim file path from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNewPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g' | sed 's/\.fna/\.fa/g')
	#Copy genome build fasta file to hisat2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNewPath"
	#Trim file extension
	buildFileNoPath=$(echo $buildFileNewPath | sed 's/\.fa//g')
	#Begin hisat2 build
	echo "Beginning hisat2 build... "
	hisat2-build -p 4 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath"
	echo "hisat2 build complete!"
else
	echo "Build folder reference_hisat2_build already exists, skipping building..."
fi

#Add run inputs to output summary file
echo "$outputsPath" > $inputOutFile
echo hisat2-build -p 4 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath" >> $inputOutFile
