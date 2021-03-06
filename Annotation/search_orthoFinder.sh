#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N search_orthoFinder_jobOutput
#$ -pe smp 8
#Script to use OrthoFinder to find orthogroups and orthologs, 
# infers rooted gene trees for all orthogroups and identifies 
# all of the gene duplication events in those gene trees
#Usage: qsub search_orthoFinder.sh analysisType proteomeFastaList
#Usage ex: qsub search_orthoFinder.sh MSA trimmed_run1E05_assemblyTrinity trimmed_run1Y05_assemblyTrinity trimmed_run1R2_assemblyTrinity trimmed_run1Y023_5_assemblyTrinity trimmed_run1PA_assemblyTrinity trimmed_run1Sierra_assemblyTrinity PA42_v4.1_proteins
#Usage ex: qsub search_orthoFinder.sh default trimmed_run1E05_assemblyTrinity/clusteredNucleotides_cdhit_0.98 trimmed_run1Y05_assemblyTrinity/clusteredNucleotides_cdhit_0.98 trimmed_run1R2_assemblyTrinity/clusteredNucleotides_cdhit_0.98 trimmed_run1Y023_5_assemblyTrinity/clusteredNucleotides_cdhit_0.98 trimmed_run1PA_assemblyTrinity/clusteredNucleotides_cdhit_0.98 trimmed_run1Sierra_assemblyTrinity/clusteredNucleotides_cdhit_0.98 PA42_v4.1_proteins
#Usage ex: qsub search_orthoFinder.sh default sortedCoordinate_samtoolsHisat2_run2E05_assemblyPA42_v3.0Trinity sortedCoordinate_samtoolsHisat2_run2Y05_assemblyPA42_v3.0Trinity sortedCoordinate_samtoolsHisat2_run2Y023_5_assemblyPA42_v3.0Trinity sortedCoordinate_samtoolsHisat2_run2R2_assemblyPA42_v3.0Trinity sortedCoordinate_samtoolsHisat2_run2PA_assemblyPA42_v3.0Trinity sortedCoordinate_samtoolsHisat2_run2Sierra_assemblyPA42_v3.0Trinity PA42_v3.0_proteins
#Usage ex: qsub search_orthoFinder.sh default sortedCoordinate_samtoolsHisat2_run2E05_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run2Y05_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run2Y023_5_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run2R2_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run2PA_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run2Sierra_assemblyPA42_v3.0Trinity/clusteredNucleotides_cdhit_0.98 PA42_v3.0_proteins
#Usage ex: qsub search_orthoFinder.sh default sortedCoordinate_samtoolsHisat2_run1E05_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run1Y05_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run1Y023_5_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run1R2_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run1PA_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 sortedCoordinate_samtoolsHisat2_run1Sierra_assemblyPA42_v4.1Trinity/clusteredNucleotides_cdhit_0.98 PA42_v4.1_proteins

#Load necessary modules
module load bio

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "No arguments supplied... exiting"
   	exit 1
fi
#Set software path
softwarePath=$(grep "orthoFinder:" ../InputData/softwarePaths.txt | tr -d " " | sed "s/orthoFinder://g")
#Set outputs absolute path
outputPath=$(grep "orthoFinder:" ../InputData/outputPaths.txt | tr -d " " | sed "s/orthoFinder://g")
dirFlag=0
runNum=1
#Make a new directory for each run
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
for i in "${@:2}"; do #Skip first two argument
	#Determine input proteome
	if [[ "$i" == *assemblyTrinity* || "$i" == *assemblyStringtie* ]]; then
		#Retrieve reads input absolute path
		inputsPath=$(grep "assemblingFree:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingFree://g")
	elif [[ "$i" == *assembly*Trinity* || "$i" == *assembly*Stringtie* ]]; then
		#Retrieve reads input absolute path
		inputsPath=$(grep "assemblingGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/assemblingGenome://g")
	elif [[ "$i" == *proteins ]]; then
		#Retrieve genome reference absolute path for querying
		inputsPath=$(grep "proteinSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/proteinSequences://g")
		#Copy input fasta files
		cp "$inputsPath" "$inputsDir"
	elif [[ "$i" == *cds ]]; then
		#Retrieve genome reference absolute path for querying
		inputsPath=$(grep "codingSequences:" ../InputData/inputPaths.txt | tr -d " " | sed "s/codingSequences://g")
		inputsPath=$(dirname "$inputsPath")
		inputsPath=$(echo "$inputsPath"/decoded_transdecoder/*transdecoder.pep)
		#Copy input fasta files
		cp "$inputsPath" "$inputsDir"
	else
		#Error message
		echo "Invalid fasta entered (species proteome expected)... exiting!"
		exit 1
	fi
	#Determine input file type
	if [[ "$i" == */clusteredNucleotide* ]]; then
		#Set inputs path
		inputsPath="$inputsPath"/"$i"/decoded_transdecoder/cdhitEst.transdecoder.pep
		#Copy input fasta files
		tag=$(echo "$i" | sed 's/\//./g')
		cp "$inputsPath" "$inputsDir"/"$tag".cdhitEst.transdecoder.pep
	elif [[ "$i" == */clusteredProtein* ]]; then
		#Set inputs path
		inputsPath="$inputsPath"/"$i"/decoded_transdecoder/cdhit.transdecoder.pep
		#Copy input fasta files
		tag=$(echo "$i" | sed 's/\//./g')
		cp "$inputsPath" "$inputsDir"/"$tag".cdhit.transdecoder.pep
	elif [[ "$i" == *proteins ]]; then
		#Copy input fasta files
		cp "$inputsPath" "$inputsDir"
	else 
		#Set inputs path
		inputsPath=$(echo "$inputsPath"/"$i"/decoded_transdecoder/*.transdecoder.pep)
		#Copy input fasta files
		cp "$inputsPath" "$inputsDir"/"$i".transdecoder.pep
	fi
done
#Move to output folder
cd "$outputFolder"
#Name output file of inputs
inputOutFile="$outputFolder"/searched_orthoFinder_summary.txt
#Use OrthoFinder to find orthologs
echo "Beginning OrthoFinder search..."
#Check OrthoFinder species tree inference method
if [[ "$1" == MSA ]]; then #Multiple sequence alignment
	"$softwarePath"/orthofinder -f "$inputsDir" -t 8 -M msa
	#Output run commands to summary file
	echo "$softwarePath"/"orthofinder -f" "$inputsDir" "-t 24 -M msa" > "$inputOutFile"
else #Default
	"$softwarePath"/orthofinder -f "$inputsDir" -t 8
	#Output run commands to summary file
	echo "$softwarePath"/"orthofinder -f" "$inputsDir" "-t 24" > "$inputOutFile"
fi
echo "OrthoFinder search complete!"
