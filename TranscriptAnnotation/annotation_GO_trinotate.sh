#!/bin/bash
#Script to use trinotate to generate an annotation report
#Usage: bash annotation_GO_trinotate.sh transcriptomeFolder
#Alternate usage Ex: bash annotation_GO_trinotate.sh PA42_cds
#Alternate usage Ex: bash annotation_GO_trinotate.sh PA42_transcripts

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "trinotatePackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/trinotatePackage://g")
if [[ "$1" == PA42_cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/codingSequencesDB://g")
elif [[ "$1" == PA42_transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/transcriptSequencesDB://g")
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Move to output folder
cd "$inputsPath"
#Set inputs path
inputsPath="$inputsPath"/"$1"
#Retrieve GO annotations
"$softsPath"/Trinotate /util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls "$inputsPath" -G --include_ancestral_terms > go_annotations.txt