#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
#Prepare for alignment
cd ..
dirFlag=0
runNum=1
genomeFile=$(head -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt")
buildFile=$(tail -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePaths.txt")
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	#Make a new directory for each alignment run
	while [ $dirFlag -eq 0 ]; do
		#Tophat output directory name
		tophatOut="aligned_tophat2_run$runNum"
		mkdir "$tophatOut"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for run $runNum of tophat2 alignment on $f1 data..."
		fi
	done
	#Name output file of inputs
	inputOutFile="$tophatOut"/"$tophatOut"_summary.txt
	#Build output directory for Tophat reference
	buildOut="reference_bowtie2_build"
	#Trim .fa file extension from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//')
	#Loop through all forward and reverse paired reads and run tophat2 on each pair
	# using 8 threads
	for f2 in "$f1"/*pForward.fq.gz; do
		#Trim extension from current file name
		curSample=$(echo $f2 | sed 's/.pForward\.fq\.gz//')
		#Trim file path from current file name
		curSampleNoPath=$(basename $f2)
		curSampleNoEx=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
		#Begin Tophat run for current sample
		echo "Sample $curSampleNoEx is being aligned..."
		tophat2 -p 8 -G "$genomeFile" -o "$tophatOut"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f2" "$curSample"_pReverse.fq.gz
		echo "Sample $curSampleNoEx has been aligned!"
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "tophat2 -p 8 -G "$genomeFile" -o "$tophatOut"/"$curSampleNoEx" "$buildOut"/"$buildFileNoEx" "$f2" "$curSample"_pReverse.fq.gz" >> $inputOutFile
	done
done
