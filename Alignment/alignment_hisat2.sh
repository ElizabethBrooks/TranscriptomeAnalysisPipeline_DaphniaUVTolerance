#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_hisat2_jobOutput
#$ -pe smp 8
#Script to perform hisat2 alignment of trimmed
# paired end reads
#Note that a hisat2 genome refernce build folder needs to be generated first
#Usage: qsub alignment_hisat2.sh alignmentTarget trimmedFolder optionalAssemblyFolder minIntronLength maxIntronLength
#Usage Ex: qsub alignment_hisat2.sh genome trimmed_run1 20 14239
#Usage Ex: qsub alignment_hisat2.sh genome trimmed_run1 20 23554
#Alternate usage Ex: qsub alignment_hisat2.sh assembly trimmed_run1 trimmed_run1E05_assemblyTrinity 20 14239
#Alternate usage Ex: qsub alignment_hisat2.sh assemblyStats trimmed_run1 trimmed_run1E05_assemblyTrinity 20 14239
#Default usage Ex: qsub alignment_hisat2.sh genome trimmed_run1

#Required modules for ND CRC servers
module load bio
#module load bio/hisat2/2.1.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == assembly* ]]; then
	analysisInput="assembly"
	#Retrieve reads input absolute path
	assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	#Retrieve build transcriptome files absolute path
	buildInputsPath="$assemblyPath"/"$3"
	#Retrieve transcriptome reference absolute path for alignment
	buildFile="$assemblyPath"/"$3"/"Trinity.fasta"
	#Retrieve alignment outputs absolute path
	outputsPath="$assemblyPath"/"$3"
	#Determine if intron lengths were entered
	if [[ -z "$4" || -z "$5" ]]; then #Arguments were not entered
		minIntron=-1
		maxIntron=-1
	else
		minIntron=$4
		maxIntron=$5
	fi
elif [[ "$1"  == genome* ]]; then
	analysisInput="trimmed"
	#Retrieve build genome files absolute path
	buildInputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/buildingGenome://g")
	#Retrieve genome reference absolute path for alignment
	buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Retrieve alignment outputs absolute path
	outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	#Determine if intron lengths were entered
	if [[ -z "$3" || -z "$4" ]]; then #Arguments were not entered
		minIntron=-1
		maxIntron=-1
	else
		minIntron=$3
		maxIntron=$4
	fi
else
	echo "ERROR: The input alignment target is not valid... exiting!"
	exit 1
fi
#Retrieve reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
trimmedFolder="$2"
#Move to outputs directory
cd "$outputsPath"
#Prepare for mapping
dirFlag=0
runNum=1
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	#Hisat output directory name
	outputFolder="aligned_hisat2_run$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of hisat2 $1 alignment on $2 data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Build output directory for Hisat reference
buildOut="$buildInputsPath"/"reference_hisat2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/')
buildFileNoEx=$(echo $buildFileNoEx | sed 's/\.fa//')
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
	if [[ $minIntron == -1 || $maxIntron == -1 ]]; then #Arguments were not entered
		#Run hisat2 with default settings
		echo "Sample $curSampleNoPath is being aligned using default settings..."
		hisat2 -p 8 -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "hisat2 -p 8 -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
	else #Run hisat2 using input intron lengths
		echo "Sample $curSampleNoPath is being aligned using input intron lengths..."
		hisat2 -p 8 --min-intronlen $minIntron --max-intronlen $maxIntron -q -x "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "hisat2 -p 8 --min-intronlen $minIntron --max-intronlen $maxIntron -q -x" "$buildOut"/"$buildFileNoEx" -1 "$f1" -2 "$curSample"_pReverse.fq.gz -S "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam --summary-file "$outputFolder"/"$curSampleNoPath"/alignedSummary.txt >> "$inputOutFile"
	fi
	echo "Sample $curSampleNoPath has been aligned!"
	#Convert or clean up bam files depending on analysis
	if [[ "$1"  == *Stats ]]; then #Clean up excess bam files, if assemly was input
		rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
	else #Convert output sam files to bam format for downstream analysis
		echo "Sample $curSampleNoPath is being converted..."
		samtools view -@ 8 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam > "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
		echo "Sample $curSampleNoPath has been converted!"
		#Remove the now converted .sam file
		rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam
		#Add run inputs to output summary file
		echo samtools view -@ 8 -bS "$outputFolder"/"$curSampleNoPath"/accepted_hits.sam ">" "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam >> "$inputOutFile"
	fi
done
#Copy previous summaries
cp "$inputsPath"/"$trimmedFolder"/*.txt "$outputFolder"
cp "$buildOut"/*.txt "$outputFolder"
