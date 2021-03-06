#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N building_bowtie2_jobOutput
#$ -q debug
#Script to generate a bowtie2 genome refernce build folder
#Usage: qsub building_bowtie2.sh trimmedOrAssemblyFolder
#Usage Ex: qsub building_bowtie2.sh trimmed_run1
#Alternate usage Ex: qsub building_bowtie2.sh trimmed_run1E05_assemblyTrinity

#Required modules for ND CRC servers
module load bio
#Prepare for alignment
dirFlag=0
runNum=1
#Determine which analysis folder was input
if [[ "$1"  == *assembly* ]]; then
	#Retrieve build outputs absolute path
	outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	outputsPath="$outputsPath"/"$1"
	#Retrieve transcriptome reference absolute path for alignment
	buildFile="$outputsPath"/"Trinity.fasta"
elif [[ "$1"  == trimmed* ]]; then
	#Retrieve genome reference absolute path for alignment
	buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Retrieve build outputs absolute path
	outputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/buildingGenome://g")
else
	echo "Input analysis target (genome or assembly folder) is not valid... exiting!"
	exit 1
fi
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
	buildFileNewPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g')
	#Copy genome build fasta file to hisat2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNewPath"
	#Trim file extension
	buildFileNoPath=$(echo $buildFileNewPath | sed 's/\.fa//g')	
	#Begin Bowtie2 build
	echo "Beginning bowtie2 build... "
	bowtie2-build "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath"
	echo "Bowtie2 build complete!"
else
	echo "Build folder reference_bowtie2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo bowtie2-build "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath" >> $inputOutFile
