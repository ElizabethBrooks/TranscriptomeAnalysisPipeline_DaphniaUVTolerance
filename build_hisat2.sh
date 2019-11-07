#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N build_bowtie2_jobOutput
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
buildOut="reference_hisat2_build"
mkdir "$buildOut"
#Name output file of inputs
inputOutFile="$buildOut"/"$buildOut"_summary.txt
if [ $? -eq 0 ]; then
	#Copy genome build fasta file to Tophat build folder
	buildFileNoPath=$(basename $buildFile)
	cp "$buildFile" "$buildOut"/"$buildFileNoPath"
	#Trim .fa file extension from build file
	buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//')
	#Begin Bowtie2 build
	echo "Beginning bowtie2 build... "
	hisat2-build -p 8 -f "$buildOut"/"$buildFileNoPath" "$buildOut"/"$buildFileNoEx"
	echo "Bowtie2 build complete!"
else
	echo "Build folder reference_bowtie2_build already exists, skipping building..."
fi
#Add run inputs to output summary file
echo "bowtie2-build "$buildOut"/"$buildFileNoPath" "$buildOut"/"$buildFileNoEx >> $inputOutFile