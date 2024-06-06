#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 4

# Script to perform hisat2 alignment of trimmed
# paired end reads
# Note that a hisat2 genome refernce build folder needs to be generated first
# usage: qsub alignment_hisat2.sh trimmedFolder
# usage Ex: qsub alignment_hisat2.sh trimmed
# usage Ex: qsub alignment_hisat2.sh trimmed_run1
# usage Ex: qsub alignment_hisat2.sh trimmed_run2

#Required modules for ND CRC servers
module load bio/2.0
#module load bio/hisat2/2.1.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve build genome files absolute path
buildInputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/buildingGenome://g")
#Retrieve genome reference absolute path for alignment
buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

#Retrieve reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
trimmedFolder=$1
#Move to outputs directory
cd "$outputsPath"

#Hisat output directory name
outputFolder="aligned_hisat2_"$2
mkdir "$outputFolder"
# check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi

#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Add software version to output summary file
hisat2 --version > $inputOutFile
#Build output directory for Hisat reference
buildOut="$buildInputsPath"/"reference_hisat2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta//' | sed 's/\.fna//' | sed 's/\.fa//')
#Loop through all forward and reverse paired reads and run Hisat2 on each pair
# using 8 threads and samtools to convert output sam files to bam
for f1 in "$inputsPath"/"$trimmedFolder"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	#Run hisat2 with default settings
	echo "Sample $curSampleNoPath is being aligned and converted..."
	hisat2 -p 4 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
	#Add run inputs to output summary file
	echo $curSampleNoPath >> $inputOutFile
	echo "hisat2 -p 4 -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
	#Convert output sam files to bam format for downstream analysis
	samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
	#Remove the now converted .sam file
	rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	#Add run inputs to output summary file
	echo samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam ">" "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam >> "$inputOutFile"
	# status message
	echo "Sample $curSampleNoPath has been aligned and converted!"
done
#Copy previous summaries
cp "$inputsPath"/"$trimmedFolder"/*.txt "$outputFolder"
cp "$buildOut"/*.txt "$outputFolder"
