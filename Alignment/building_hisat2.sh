#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N building_hisat2_jobOutput
#$ -pe smp 8
#$ -q debug
#Script to generate a hisat2 genome refernce build folder
#Usage: qsub building_hisat2.sh trimmedOrAssemblyFolder
#Usage ex: qsub building_hisat2.sh trimmed_run1
#Alternate usage ex: qsub building_hisat2.sh trimmed_run1E05_assemblyTrinity
#Alternate usage ex: qsub building_hisat2.sh sortedCoordinate_samtoolsHisat2_run2E05_assemblyPA42_v3.0Trinity
#Alternate usage ex: qsub building_hisat2.sh sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98

#Required modules for ND CRC servers
module load bio
#module load bio/hisat2/2.1.0
#Prepare for alignment
dirFlag=0
runNum=1
#Determine which analysis folder was input
#Determine the type of assembly
if [[ "$1" == *assemblyTrinity* || "$1" == *assemblyStringtie* ]]; then
	#Retrieve reads input absolute path
	outputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	#Retrieve build outputs absolute path
	outputsPath="$outputsPath"/"$1"
	#Retrieve transcriptome reference absolute path for alignment
	if [[ "$1" == *clusteredNucleotide* ]]; then
		buildFile="$outputsPath"/"cdhitEst"
	else
		buildFile="$outputsPath"/"Trinity.fasta"
	fi
elif [[ "$1" == *assembly*Trinity* || "$1" == *assembly*Stringtie* ]]; then
	#Retrieve reads input absolute path
	outputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	#Retrieve build outputs absolute path
	outputsPath="$outputsPath"/"$1"
	#Retrieve transcriptome reference absolute path for alignment
	if [[ "$1" == *clusteredNucleotide* ]]; then
		buildFile="$outputsPath"/"cdhitEst"
	else
		buildFile="$outputsPath"/"Trinity.fasta"
	fi
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
	buildFileNewPath=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/g' | sed 's/\.fna/\.fa/g')
	#Copy genome build fasta file to hisat2 build folder
	cp "$buildFile" "$outputFolder"/"$buildFileNewPath"
	#Trim file extension
	buildFileNoPath=$(echo $buildFileNewPath | sed 's/\.fa//g')
	#Begin hisat2 build
	echo "Beginning hisat2 build... "
	hisat2-build -p 8 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath"
	echo "hisat2 build complete!"
else
	echo "Build folder reference_hisat2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo "$outputsPath" > $inputOutFile
echo hisat2-build -p 8 -f "$outputFolder"/"$buildFileNewPath" "$outputFolder"/"$buildFileNoPath" >> $inputOutFile
