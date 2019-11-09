#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N build_bowtie2_jobOutput
#Required modules for ND CRC servers
module load bio
#Prepare for alignment
cd ..
dirFlag=0
runNum=1
buildFile=$(tail -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt")
#Build reference genome if folder does not exist
#Build output directory
buildOut="reference_bowtie2_build"
mkdir "$buildOut"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	#Build files already exsist
	echo "Build files for bowtie2 already found in $buildOut... exiting"
	exit 1
else
	#Build files do not exsist
	echo "Creating $buildOut folder for bowtie2 build..."
fi
#Name output file of inputs
inputOutFile="$buildOut"/"$buildOut"_summary.txt
if [ $? -eq 0 ]; then
	#Trim file path from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNoPath=$(echo $buildFileNoPath | sed 's/\.fasta/.fa/g')
	#Copy genome build fasta file to bowtie2 build folder
	cp "$buildFile" "$buildOut"/"$buildFileNoPath"
	#Begin Bowtie2 build
	echo "Beginning bowtie2 build... "
	bowtie2-build "$buildOut"/"$buildFileNoPath" "$buildOut"/"$buildFileNoPath"
	echo "Bowtie2 build complete!"
else
	echo "Build folder reference_bowtie2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo bowtie2-build "$buildOut"/"$buildFileNoPath" "$buildOut"/"$buildFileNoPath" >> $inputOutFile