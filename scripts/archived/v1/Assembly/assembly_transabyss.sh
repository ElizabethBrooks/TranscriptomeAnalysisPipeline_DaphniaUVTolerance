#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_transabyss_jobOutput
#$ -pe smp 8
#Script to perform de novo transcriptome assembly using Transabyss
#Usage: qsub assembly_transabyss.sh trimmedFolder genotype
#Usage Ex: qsub assembly_transabyss.sh trimmed_run1 E05
#Note that the genome version input is for output file naming purposes only

#Required modules for ND CRC servers
module load bio/transabyss
#module load python
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
if [[ "$1"  != trimmed* ]]; then
	echo "ERROR: The $1 folder of aligned bam files were not found... exiting"
	exit 1
fi
#Retrieve trimmed reads input absolute path
inputsPath=$(grep "trimming:" ../InputData/outputPaths.txt | tr -d " " | sed "s/trimming://g")
#Retrieve outputs path
outputsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
#Move to outputs directory
outputFolder="$outputsPath"/"$1""$2"_assemblyTrinity
mkdir "$outputFolder"
#Re-set reads file paths using input genotype tag
sed "s/GENEOTYPE/$2/g" "$samplesPath" > "$outputFolder"/tmpSamplesFile.txt
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
cd "$outputsPath"
#TO DO: Need igraph (pip install python-igraph)
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Loop through all forward and reverse paired reads and
# store the file names in an array list
#Set the flag for paired-end sample input to transabyss
sampleList="--pe"
let counter=$counter+1
for f1 in "$inputsPath"/"$1"/*"$2"*pForward.fq.gz; do
	#Store current sample file name
	sampleList=$sampleList" "$f1
done
#Begin transabyss assembly
transabyss --threads 8 $sampleList --SS --name "$2" --outdir "$outputFolder"
echo "Transabyss assembly of $1 data is complete!"
#Add run inputs to output summary file
echo "$curSampleNoPath" > "$inputOutFile"
echo "transabyss --threads 8 $sampleList --SS --name $2 --outdir $outputFolder" > "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"