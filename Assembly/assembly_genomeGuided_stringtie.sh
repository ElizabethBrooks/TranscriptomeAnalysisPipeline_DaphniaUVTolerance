#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assemblyGenome_stringtie_jobOutput
#$ -pe smp 8
#Script to assemble transcripts using a reference genome and stringtie
#Usage: qsub assembly_genomeGuided_stringtie.sh sortedFolder genotype genome
#Usage Ex: qsub assembly_genomeGuided_stringtie.sh sortedCoordinate_samtoolsHisat2_run1 E05 PA42_v4.1

##Retrieve software absolute path
softwarePath=$(grep "stringtie:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/stringtie://g")
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine if the folder name was input in the correct format
if [[ "$1" == *\/* ]] || [[ "$1" == *\\* ]]; then
	echo "ERROR: Please enter folder names without a trailing forward slash (/)... exiting"
	exit 1
fi
#Determine if the correct analysis folder was input
if [[ "$1"  != sortedCoordinate* ]]; then
	echo "ERROR: The "$1" folder of coordinate sorted aligned bam files were not found... exiting"
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
#Retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#Retrieve genome reference absolute path for alignment
genomeFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")
#Retrieve alignment outputs absolute path
outputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
#Create output directory
outputFolder="$outputsPath"/"$1""$2"_assembly"$3"Stringtie
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to outputs directory
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/"$1""$2"_assembly"$3"Stringtie_summary.txt
sampleList=""
#Run stringtie on each coordinate-sorted bam files using 8 threads, and a maximum intron
# length that makes most sense given your targeted organism
for f1 in "$inputsPath"/"$1"/*_"$2"_*/*.bam; do
	#Retrieve curent sample
	curAlignedSample="$f1"
	#Trim file path from current folder name
	curSampleNoPath=$(echo $f1 | sed 's/accepted\_hits\.bam//g')
	curSampleNoPath=$(basename "$curSampleNoPath")
	#The main input of the program is a BAM file with RNA-Seq 
	# read mappings which must be sorted by their genomic location 
	echo "Beginning assembly of $curSampleNoPath..."
	"$softwarePath"/stringtie "$curAlignedSample" -p 8 -G "$genomeFile" -C "$curSampleNoPath".cov_refs.gtf -o "$curSampleNoPath".stringtie.gtf
	echo "Assembly of $curSampleNoPath complete!"
	#Add sample to list
	sampleList="$curSampleNoPath".stringtie.gtf "$sampleList"
	#Add run inputs to output summary file
	echo "stringtie "$curAlignedSample" -p 8 -G "$genomeFile" -C "$curSampleNoPath".cov_refs.gtf -o "$curSampleNoPath".stringtie.gtf" > "$inputOutFile"
done
#Merge to generate a non-redundant set of transcripts observed in any of the reads
"$softwarePath"/stringtie "$sampleList" -p 8 -G "$genomeFile" -o "$1""$2"_assembly"$3"Stringtie.gtf --merge 
#Add run inputs to output summary file
echo "stringtie "$sampleList" -p 8 -G "$genomeFile" -o "$1""$2"_assembly"$3"Stringtie.gtf --merge" > "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"