#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_genomeGuided_trinity_jobOutput
#$ -pe smp 8
#Script to perform genome-guided Trinity de novo transcriptome assembly
#Usage: qsub assembly_genomeGuided_trinity.sh sortedFolder genotype maxIntronLength
#Usage Ex: qsub assembly_genomeGuided_trinity.sh sortedCoordinate_samtoolsTophat2_run1 Sierra 14239

#Required modules for ND CRC servers
module load bio/2.0
#Retrieve aligned reads input absolute path
inputsPath=$(grep "sorting:" ../InputData/outputPaths.txt | tr -d " " | sed "s/sorting://g")
#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
#Retrieve assembly outputs absolute path
outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Create output directory
outputFolder="$outputsPath"/"$1""$2"_assembly_Trinity
mkdir "$outputFolder"
#Move to outputs directory
cd "$outputFolder"
#Prepare for analysis
dirFlag=0
runNum=1
COUNTER=0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != sortedCoordinate* ]]; then
	echo "ERROR: The "$1" folder of aligned bam files were not found... exiting"
	exit 1
fi
#Determine what analysis method was used for the input folder of data
if [[ "$1" == *"Hisat2"*  ]]; then
	#Set analysis method for folder naming
	analysisMethod="Hisat2"
elif [[ "$1" == *"Tophat2"* ]]; then
	#Set analysis method for folder naming
	analysisMethod="Tophat2"
else
	echo "ERROR: The sorted "$1" folder of bam files were not found... exiting"
	exit 1
fi
#Name output file of inputs
inputOutFile="$outputFolder"/"$1""$2"_assembly_summary.txt
#Merge and re-coordinate sort the set of bam files
readFiles=$(echo "$inputsPath"/"$1"/*"$2"*/*.bam)
mergedBam="$outputFolder"/"$1""$2"_merged.bam
samtools merge -@ 8 "$mergedBam" "$readFiles"
sortedBam="$outputFolder"/"$1""$2"_mergedSorted.bam
samtools sort -@ 8 -o "$sortedBam" "$mergedBam"
rm "$mergedBam"
#Run Trinity on coordinate-sorted bam files using 8 threads, and a maximum intron
# length that makes most sense given your targeted organism
echo "Sample $curSampleNoPath is being assembled..."
Trinity --genome_guided_bam "$sortedBam" --genome_guided_max_intron "$3" --max_memory 50G --CPU 8
rm "$sortedBam"
#Add run inputs to output summary file
echo "samtools merge --threads 8" "$mergedBam" "$readFiles" > "$inputOutFile"
echo "samtools sort -@ 8 -o" "$sortedBam" "$mergedBam" >> "$inputOutFile"
echo "Trinity --genome_guided_bam" "$f1" "--genome_guided_max_intron" "$3" "--max_memory 50G --CPU 8" >> "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"