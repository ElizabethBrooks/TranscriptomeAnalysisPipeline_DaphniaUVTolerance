#!/bin/bash
#Script to generate a multi FASTA file for all transcripts in a GFF file
#Usage: bash generateTranscripts_cufflinks.sh sortedNameFolder analysisTarget
#Usage Ex: bash generateTranscripts_cufflinks.sh sortedCoordinate_samtoolsHisat2_run3 variantCallingBcftools_filteredMapQ
#Usage Ex: bash generateTranscripts_cufflinks.sh sortedCoordinate_samtoolsHisat2_run3 variantCalling_filteredMapQ
#Usage Ex: bash generateTranscripts_cufflinks.sh genomeReference

#Required modules for ND CRC servers
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve genome reference and features absolute paths
genomeFeatFile=$(grep "genomeFeatures:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeFeatures://g")

#Determine what analysis method was used for the input folder of data
if [[ "$1" == sorted* ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
	type=$(echo "$2" | cut -d"_" -f2)
	inputsPath="$inputsPath"/"$1"/"$2"
	inputFeatFile="$inputsPath"/"$type"_consensusFeatures.gff
	inputsPath="$inputsPath"/"$type"_consensus.fa
elif [[ "$1" == genomeReference ]]; then
	#Retrieve sorted reads input absolute path
	inputsPath=$(grep "genomeReference:" ../InputData/inputPaths.txt | tr -d " " | sed "s/genomeReference://g")
else
	echo "ERROR: Invalid sorted folder of bam files entered... exiting"
	exit 1
fi

#Name output file of inputs
outputsPath=$(dirname "$inputsPath")
inputOutFile="$outputsPath"/generateTranscipts_summary.txt
errorOut="$outputsPath"/generateTransciptsError_summary.txt

#Generate a fasta index file using samtools
samtools faidx "$inputsPath"

#Generate a FASTA file with the DNA sequences for all transcripts in the GFF file
gffread -w "$outputsPath"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile" 2> "$errorOut"
echo "Generate fasta file with the DNA sequences for all transcripts in the updated GFF file" > "$inputOutFile"
echo gffread -w "$outputsPath"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile" >> "$inputOutFile"

#Find each error coordinate and replace with the correct one
if [[ "$1" == sorted* ]]; then
	while read -r line; do
		find=$(cut -d" " -f5 "$line" | sed "s/)//g" | sed "s/(//g")
		replace=$(cut -d" " -f12 "$line")
		#Output status message
		echo "Find: $find & Replace: $replace"
		#Update coordinates in feature file
		sed -i "s/\t$find\t/\t$replace\t/g" "$inputFeatFile"
		#Generate a FASTA file with the DNA sequences for all transcripts in the updated GFF file
		gffread -w "$outputsPath"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile"
		echo "Generate fasta file with the DNA sequences for all transcripts in the updated GFF file" >> "$inputOutFile"
		echo gffread -w "$outputsPath"/transcripts_cufflinks.fa -g "$inputsPath" "$genomeFeatFile" >> "$inputOutFile"
	done < "$errorOut"
fi