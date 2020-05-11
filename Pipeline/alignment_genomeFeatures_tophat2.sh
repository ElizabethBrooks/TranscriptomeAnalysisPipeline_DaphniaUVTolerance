#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N alignment_tophat2_jobOutput
#$ -pe smp 8
#Script to perform tophat2 alignment of trimmed
# paired end reads
#Note that a bowtie2 genome refernce build folder needs to be generated first
#Usage: qsub alignment_tophat2.sh alignmentTarger trimmedFolder optionalAssemblyFolder minIntronLength maxIntronLength
#Usage Ex: qsub alignment_tophat2.sh genome trimmed_run1 4 14239
#Alternate usage Ex: qsub alignment_tophat2.sh assembly trimmed_run1E05_assemblyTrinity 20 14239
#Default usage Ex: qsub alignment_tophat2.sh genome trimmed_run1

#Required modules for ND CRC servers
module load bio
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == assembly ]]; then
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
		minIntron=$3
		maxIntron=$4
	else
		minIntron=-1
		maxIntron=-1
	fi
elif [[ "$1"  == genome ]]; then
	analysisInput="trimmed"
	#Retrieve build genome files absolute path
	buildInputsPath=$(grep "buildingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/building://g")
	#Retrieve genome reference absolute path for alignment
	buildFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
	#Retrieve alignment outputs absolute path
	outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligning://g")
	#Determine if intron lengths were entered
	if [[ -z "$3" || -z "$4" ]]; then #Arguments were not entered
		minIntron=$3
		maxIntron=$4
	else
		minIntron=-1
		maxIntron=-1
	fi
else
	echo "ERROR: The input folder of trimmed or assembled files were not found... exiting"
	exit 1
fi
#Retrieve reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
#Retrieve genome features absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Move to outputs directory
cd "$outputsPath"
#Generate tronscriptome index files
tophat -G "$genomeFile" --transcriptome-index="$outputsPath"/transcriptome_data/known pa42
#Prepare for alignment
dirFlag=0
runNum=1
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	#Tophat output directory name
	outputFolder="aligned_tophat2_run$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of tophat2 alignment on $1 data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
trimmedFolder="$2"
#Build output directory for Tophat reference
buildOut="$buildInputsPath"/"reference_bowtie2_build"
#Trim .fa file extension from build file
buildFileNoPath=$(basename $buildFile)
buildFileNoEx=$(echo $buildFileNoPath | sed 's/\.fasta/\.fa/')
buildFileNoEx=$(echo $buildFileNoEx | sed 's/\.fa//')
#Loop through all forward and reverse paired reads and run tophat2 on each pair
# using 8 threads
for f1 in "$inputsPath"/"$trimmedFolder"/*pForward.fq.gz; do
	#Trim extension from current file name
	curSample=$(echo $f1 | sed 's/.pForward\.fq\.gz//')
	#Trim file path from current file name
	curSampleNoPath=$(basename $f1)
	curSampleNoPath=$(echo $curSampleNoPath | sed 's/.pForward\.fq\.gz//')
	#Determine if intron lengths were entered
	if [[ $minIntron == -1 || $maxIntron == -1 ]]; then #Arguments were not entered
		#Run tophat2 with default settings
		echo "Sample $curSampleNoPath is being aligned..."
		tophat2 -p 8 --transcriptome-index="$outputsPath"/transcriptome_data/known pa42 -o "$outputFolder"/"$curSampleNoPath" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "tophat2 -p 8 --transcriptome-index=$outputsPath/transcriptome_data/known pa42 -o $outputFolder/$curSampleNoPath $buildOut/$buildFileNoEx $f1 $curSample_pReverse.fq.gz" >> "$inputOutFile"
	else #Run tophat2 using input intron lengths
		echo "Sample $curSampleNoPath is being aligned..."
		tophat2 -p 8 --transcriptome-index="$outputsPath"/transcriptome_data/known pa42 -i $minIntron -I $maxIntron -G "$genomeFile" -o "$outputFolder"/"$curSampleNoPath" "$buildOut"/"$buildFileNoEx" "$f1" "$curSample"_pReverse.fq.gz
		#Add run inputs to output summary file
		echo $curSampleNoPath >> $inputOutFile
		echo "tophat2 -p 8 --transcriptome-index=$outputsPath/transcriptome_data/known pa42 -i $minIntron -I $maxIntron -G $genomeFile -o $outputFolder/$curSampleNoPath $buildOut/$buildFileNoEx $f1 $curSample_pReverse.fq.gz" >> "$inputOutFile"
	fi
	echo "Sample $curSampleNoPath has been aligned!"
	#Clean up excess alignment files, if assemly was input
	if [[ "$1"  == assembly ]]; then
		rm "$outputFolder"/"$curSampleNoPath"/accepted_hits.bam
		rm "$outputFolder"/"$curSampleNoPath"/deletions.bed
		rm "$outputFolder"/"$curSampleNoPath"/junctions.bed
		rm "$outputFolder"/"$curSampleNoPath"/prep_reads.info
		rm "$outputFolder"/"$curSampleNoPath"/insertions.bed
		rm -r "$outputFolder"/"$curSampleNoPath"/logs
		rm "$outputFolder"/"$curSampleNoPath"/unmapped.bam
	fi
done
#Copy previous summaries
cp "$inputsPath"/"$trimmedFolder"/*.txt "$outputFolder"
cp "$buildOut"/*.txt "$outputFolder"
