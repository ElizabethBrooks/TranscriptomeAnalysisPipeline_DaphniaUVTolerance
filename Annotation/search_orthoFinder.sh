#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_orthoFinder_jobOutput
#$ -pe smp 24
#Script to use OrthoFinder to find orthogroups and orthologs, 
# infers rooted gene trees for all orthogroups and identifies 
# all of the gene duplication events in those gene trees
#Usage: qsub search_orthoFinder.sh proteomeFastaList
#Usage Ex: qsub search_orthoFinder.sh trimmed_run1E05_assemblyTrinity trimmed_run1Y05_assemblyTrinity trimmed_run1R2_assemblyTrinity trimmed_run1Y023_5_assemblyTrinity trimmed_run1PA_assemblyTrinity trimmed_run1Sierra_assemblyTrinity PA42_proteins

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No folder name(s) supplied... exiting"
   	exit 1
fi
#Set software path
softwarePath=$(grep "orthoFinder:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/orthoFinder://g")
#Set outputs absolute path
outputPath=$(grep "orthoFinder:" ../InputData/outputPaths.txt | tr -d " " | sed "s/orthoFinder://g")
dirFlag=0
runNum=1
#Make a new directory for each alignment run
while [ $dirFlag -eq 0 ]; do
	#Hisat output directory name
	outputFolder="$outputPath"/"searched_orthoFinder_run$runNum"
	mkdir "$outputFolder"
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of orthoFinder searching..."
	fi
done
#Loop through all input proteomes and build directory of inputs
inputsDir="$outputFolder"
for i in "$@"; do
	#Determine input proteome
	if [[ "$i" == *assembly* ]]; then
		#Retrieve reads input absolute path
		assemblyPath=$(grep "assembling:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assembling://g")
		inputsPath="$assemblyPath"/"$i"/decoded_transdecoder/Trinity.fasta.transdecoder.pep
		cp "$inputsPath" "$inputsDir"/"$i"_Trinity.fasta.transdecoder.pep
	elif [[ "$i" == PA42_proteins ]]; then
		#Retrieve genome reference absolute path for querying
		inputsPath=$(grep "proteinSequencesDB:" ../InputData/databasePaths.txt | tr -d " " | sed "s/proteinSequencesDB://g")
		cp "$inputsPath" "$inputsDir"
	else
		#Error message
		echo "Invalid fasta entered (species proteome expected)... exiting!"
		exit 1
	fi
done
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/searched_orthoFinder_summary.txt
#Use OrthoFinder to find orthologs
echo "Beginning OrthoFinder search..."
"$softwarePath"/orthofinder -f "$inputsDir" -t 24
echo "OrthoFinder search complete!"
#Output run commands to summary file
echo "$softwarePath"/"orthofinder -f" "$inputsDir" "-t 24" > "$inputOutFile"