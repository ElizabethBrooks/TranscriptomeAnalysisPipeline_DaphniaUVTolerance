#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N annotation_sanspanz_jobOutput
#Script to predict coding regions from a transcript fasta file
# using Transdecoder with genome reference and features files
#Usage: qsub annotation_sanspanz.sh decodedFolder optionalSpecies
#Usage Ex: qsub annotation_sanspanz.sh aligned_hisat2_run1/decoded_transdecoder_genomeBased
#Usage Ex: qsub annotation_sanspanz.sh trimmed_run1E05_assemblyTrinity/decoded_transdecoder E05
#Usage Ex: qsub annotation_sanspanz.sh sortedCoordinate_samtoolsTophat2_run1E05_assemblyGenomeTrinity/decoded_transdecoder E05

#Load necessary modules for ND CRC servers
module load bio/2.0
#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Determine which analysis folder was input
if [[ "$1"  == *assembly* ]]; then
	#Retrieve assembly absolute path
	inOutsPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
	species="Daphnia $2"
elif [[ "$1"  == aligned* ]]; then
	#Retrieve alignment absolute path
	inOutsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	species="Daphnia"
else
	echo "ERROR: The input alignment target is not valid... exiting!"
	exit 1
fi
#Retrieve sanspanz software package path
softsPath=$(grep "packageSanspanz:" ../InputData/inputPaths.txt | tr -d " " | sed "s/packageSanspanz://g")
#Set outputs absolute path
inOutsPath="$inOutsPath"/"$1"
#Set outputs folder
outputFolder="$inOutsPath"/annotated_sanspanz
mkdir "$outputFolder"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputFolder directory already exsists... please remove before proceeding."
	exit 1
fi
#Move to output folder
cd "$outputFolder"
#Run sanspanz for functional annotation
python "$softsPath"/runsanspanz.py -R -m SANS -s "$species" -i "$inOutsPath"/transcripts.fasta.transdecoder.pep -o sans.tab
python "$softsPath"/runsanspanz.py -R -m Pannzer -i sans.tab -f tab -o ",DE.out1,GO.out1,anno.out1"