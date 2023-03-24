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

#Required modules for ND CRC servers
module load bio
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
trimmedFolder="$1"
#Move to outputs directory
cd "$outputsPath"

#Hisat output directory name
outputFolder="aligned_hisat2"
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
	#Determine if intron lengths were entered
	if [[ "$maxIntron" == -1 ]]; then #Arguments were not entered
		if [[ "$dtAnalysis" == -1 ]]; then #Arguments were not entered
			#Run hisat2 with default settings
			echo "Sample $curSampleNoPath is being aligned using default settings..."
			hisat2 -p 4 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
			#Add run inputs to output summary file
			echo $curSampleNoPath >> $inputOutFile
			echo "hisat2 -p 4 -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
		else
			#Run hisat2 with default settings
			echo "Sample $curSampleNoPath is being aligned using default settings and dta..."
			hisat2 -p 4 --dta -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
			#Add run inputs to output summary file
			echo $curSampleNoPath >> $inputOutFile
			echo "hisat2 -p 4 --dta -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
		fi
	else #Run hisat2 using input intron lengths
		if [[ "$dtAnalysis" == -1 ]]; then #Arguments were not entered
			echo "Sample $curSampleNoPath is being aligned using input max intron length..."
			hisat2 -p 4 --max-intronlen $maxIntron -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
			#Add run inputs to output summary file
			echo $curSampleNoPath >> $inputOutFile
			echo "hisat2 -p 4 --max-intronlen $maxIntron -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
		else
			echo "Sample $curSampleNoPath is being aligned using input max intron length and dta..."
			hisat2 -p 4 --dta --max-intronlen $maxIntron -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
			#Add run inputs to output summary file
			echo $curSampleNoPath >> $inputOutFile
			echo "hisat2 -p 4 --dta --max-intronlen $maxIntron -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
		fi
	fi
	echo "Sample $curSampleNoPath has been aligned!"
	#Convert or clean up bam files depending on analysis
	if [[ "$1"  == *Stats ]]; then #Clean up excess bam files, if assemly was input
		rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	else #Convert output sam files to bam format for downstream analysis
		echo "Sample $curSampleNoPath is being converted..."
		samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
		echo "Sample $curSampleNoPath has been converted!"
		#Remove the now converted .sam file
		rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
		#Add run inputs to output summary file
		echo samtools view -@ 4 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam ">" "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam >> "$inputOutFile"
	fi
done
#Copy previous summaries
cp "$inputsPath"/"$trimmedFolder"/*.txt "$outputFolder"
cp "$buildOut"/*.txt "$outputFolder"
