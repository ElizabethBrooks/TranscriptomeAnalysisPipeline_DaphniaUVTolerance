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
runNum=0
genomeFile=$(head -n 1 "TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/genomeFilePath.txt")
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Build reference genome if folder does not exist
mkdir aligned_tophat2_build
if [ $? -eq 0 ]; then
	echo "Beginning bowtie2 build... "
	bowtie2-build /afs/crc.nd.edu/group/hoth/echo_base/genome/Daphnia_pulex.allmasked.fa aligned_tophat2_build/Daphnia_pulex.allmasked
	echo "Bowtie2 build complete!"
else
	echo "Build already exists, skipping building..."
fi
#Retrieve folders to analyze from the input arguments
for f1 in "$@"; do
	#Make a new directory for each alignment run
	while [ $dirFlag -eq 0 ]; do
		mkdir aligned_tophat2_run"$runNum"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for $runNum run of tophat2 alignment on $f1 data..."
		fi
	done
	#Loop through all forward and reverse paired reads and run tophat2 on each pair
	# using 8 threads
	mkdir aligned_tophat2_run"$runNum"/out
	for f2 in "$f1"/*pForward.fq.gz; do
		#Trim extension from current file name
		curFile=$(echo $f2 | sed 's/pForward\.fq\.gz//')
		#Trim file path from current file name
		curFileNoPath=$(basename $f2)
		curFileNoPath=$(echo $curFileNoPath | sed 's/pForward\.fq\.gz//')
		echo "Sample $curFileNoPath is being aligned..."
		tophat2 -p 8 -G "$genomeFile" -o aligned_tophat2_run"$runNum"/out aligned_tophat2_build/Daphnia_pulex.allmasked "$f2" "$curSample"pReverse.fq.gz
		echo "Sample $curFileNoPath has been aligned!"
	done
done
