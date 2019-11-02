#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N build_bowtie2_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
#Prepare for alignment
cd ..
dirFlag=0
runNum=0
buildFile=$(tail -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt")
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Build reference genome if folder does not exist
mkdir reference_bowtie2_build
if [ $? -eq 0 ]; then
	#Build output directory
	buildOut="reference_bowtie2_build"
	#Copy genome build fasta file to Tophat build folder
	buildFileNoPath=$(basename $buildFile)
	cp "$buildFile" "$buildOut"/"$buildFileNoPath"
	#Trim .fa file extension from build file
	buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fa//')
	#Begin Bowtie2 build
	echo "Beginning bowtie2 build... "
	bowtie2-build "$buildOut"/"$buildFileNoPath" "$buildOut"/"$buildFileNoEx"
	echo "Bowtie2 build complete!"
else
	echo "Build folder reference_bowtie2_build already exists, skipping building..."
fi
