#!/bin/bash
#Script to use trinotate to generate an annotation report
#Usage: bash annotation_GO_trinotate.sh transcriptomeFolder
#Alternate usage Ex: bash annotation_GO_trinotate.sh PA42_v4.1_cds
#Alternate usage Ex: bash annotation_GO_trinotate.sh PA42_v3.0_transcripts

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve Trinotate software path
softsPath=$(grep "trinotatePackage:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/trinotatePackage://g")
if [[ "$1" == *cds ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
elif [[ "$1" == *transcripts ]]; then
	#Retrieve genome reference absolute path for querying
	inputsPath=$(grep "transcriptSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/transcriptSequences://g")
else
	#Error message
	echo "Invalid fasta entered (assembled transcriptome expected)... exiting!"
	exit 1
fi
#Move to output folder
inputsPath=$(dirname $inputsPath)
cd "$inputsPath"/annotated_trinotate
#Set inputs path
inputsPath="$inputsPath"/annotated_trinotate/trinotate_annotation_report.xls
#Retrieve GO annotations
"$softsPath"/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls "$inputsPath" -G --include_ancestral_terms > go_annotations.txt