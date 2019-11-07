#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 8
#Required modules for ND CRC servers
module load bio
module load bio/hisat2/2.1.0
#Prepare for mapping
cd ..
dirFlag=0
runNum=1
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
		#Hisat output directory name
		hisatOut="aligned_hisat2_run$runNum"
		mkdir "$hisatOut"
		#Check if the folder already exists
		if [ $? -ne 0 ]; then
			#Increment the folder name
			let runNum+=1
		else
			#Indicate that the folder was successfully made
			dirFlag=1
			echo "Creating folder for run $runNum of hisat2 alignment on $f1 data..."
		fi
	done
	#Name output file of inputs
	inputOutFile="$hisatOut"/"$hisatOut"_summary.txt
	#Build output directory for Hisat reference
	buildOut="reference_hisat2_build"
	#Trim .fa file extension from build file
	buildFileNoPath=$(basename $buildFile)
	buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//')
	#Loop through all forward and reverse paired reads and run Hisat2 on each pair
	# using 8 threads and samtools to convert output sam files to bam
	for f2 in "$f1"/*pForward.fq.gz; do
		#Trim extension from current file name
		curSample=$(echo $f2 | sed 's/.pForward\.fq\.gz//')
		#Trim file path from current file name
		curSampleNoPath=$(basename $f2)
		curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
		echo "Sample $curSampleNoPath is being aligned..."
		#Create directory for current sample outputs
		mkdir "$hisatOut"/"$curSampleNoPath"
		hisat2 -p 8 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f2" -2 "$curSample"_pReverse.fq.gz -S "$hisatOut"/"$curSampleNoPath"/"$curSampleNoPath".sam --summary-file "$hisatOut"/"$curSampleNoPath"/alignedSummary.txt
		#Convert output sam files to bam format for downstream analysis
		echo "Sample $curSampleNoPath is being converted..."
		samtools view -@ 8 -bS "$hisatOut"/"$curSampleNoPath"/"$curSampleNoPath".sam > "$hisatOut"/"$curSampleNoPath"/"$curSampleNoPath".bam
		echo "Sample $curSampleNoPath has been aligned and converted!"
		#Remove the now converted .sam file
		rm "$hisatOut"/"$curSampleNoPath"/"$curSampleNoPath".sam
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "hisat2 -p 8 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f2" -2 "$curSample"_pReverse.fq.gz -S "$hisatOut"/"$curSampleNoPath"/"$curSampleNoPath".sam --summary-file "$hisatOut"/"$curSampleNoPath"/alignedSummary.txt" >> $inputOutFile
	done
done
