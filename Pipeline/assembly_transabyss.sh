#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N assembly_transabyss_jobOutput
#$ -pe smp 8
#Script to perform de novo transcriptome assembly using Transabyss
#Usage: qsub assembly_transabyss.sh trimmedFolder genomeVersion
#Usage Ex: qsub assembly_transabyss.sh trimmed_run1 PA42_v3.0
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
cd "$outputsPath"
#TO DO: Need igraph (pip install python-igraph)
#Prepare for alignment
dirFlag=0
runNum=1
counter=0
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	#Tophat output directory name
	outputFolder="assembled_transabyss_run$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of transabyss assembly on $1 data..."
	fi
done
#Name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
#Loop through all forward and reverse paired reads and
# store the file names in an array list
#Set the flag for paired-end sample input to transabyss
SARRAY[counter]="--pe"
let counter=$counter+1
for f1 in "$inputsPath"/"$1"/*pForward.fq.gz; do
	#Store current sample file name in a temp txt file
	SARRAY[counter]=" $f1"
	#Incrememnt counter
	let counter=$counter+1
done
#Begin transabyss assembly
transabyss --threads 8 $sampleList --SS --name "$2" --outdir "$outputFolder"
echo "Transabyss assembly of $1 data is complete!"
#Add run inputs to output summary file
echo "$curSampleNoPath" > "$inputOutFile"
echo "transabyss --threads 8 --pe $sampleList --SS --name $2 --outdir $outputFolder" > "$inputOutFile"
#Copy previous summaries
cp "$inputsPath"/"$1"/*.txt "$outputFolder"