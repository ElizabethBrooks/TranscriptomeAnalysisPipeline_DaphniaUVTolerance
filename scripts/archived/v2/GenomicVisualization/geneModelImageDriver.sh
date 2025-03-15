#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N geneModelImages_jobOutput
#Script to generate images for each single gene model
#Usage: qsub geneModelImageDriver.sh sortedNameFolder bamType
#Usage Ex: qsub geneModelImageDriver.sh sortedCoordinate_samtoolsHisat2_run3 filteredMapQ
#Usage Ex: qsub geneModelImageDriver.sh sortedCoordinate_samtoolsHisat2_run3 accepted_hits

#Load module for servers
module load bio

#Retrieve genome features absolute path for alignment
genomeFile=$(grep "cleanedFeatures" ../InputData/inputPaths.txt | tr -d " " | sed "s/cleanedFeatures://g")
#Retrieve outputs path
outDir=$(grep "geneModels" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneModels://g")
outDir="$outDir"_"$2"
#Make output directory
mkdir "$outDir"

#Retrieve aligned files path
bamPath=$(grep "aligningGenome" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
bamPath="$bamPath"/"$1"

#Set input bam file paths
inputBam1="$bamPath"/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/"$2".bam
inputBam2="$bamPath"/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_UV/"$2".bam
inputBam3="$bamPath"/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_UV/"$2".bam
inputBam4="$bamPath"/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_VIS/"$2".bam
inputBam5="$bamPath"/140327_I481_FCC3P1PACXX_L3_Pool_2_PA_VIS/"$2".bam
inputBam6="$bamPath"/140327_I481_FCC3P1PACXX_L4_Pool_3_PA_VIS/"$2".bam

#Set output file header
outIDs="$outDir"/multipleModels_geneIDs.csv
echo "geneID" > "$outIDs"

#Retrieve all gene IDs
grep "Name=" "$genomeFile" | cut -f 9 | cut -d "=" -f 3 > "$outDir"/tmpIDs.csv

#Loop over each gene ID
while IFS=, read -r geneName
do
	#Output status message
	echo "Processing gene $geneName..."

	#Count the number of lines for each gene and the first mRNA
	numRNA=$(($(grep -w "$geneName" "$genomeFile" | grep "mRNA-1" | wc -l)+1))
	numGene=$(grep -w "$geneName" "$genomeFile" | wc -l)

	#Check if there is more than one predicted gene model
	if [[ $numRNA == $numGene ]]; then
		#Set output gene gff path
		geneGff="$outDir"/"$geneName"_geneModel.gff 
		#Retrieve assocaited gff lines
		grep -w "$geneName" "$genomeFile" > "$geneGff"
		#Retrieve input values
		scaffoldName=$(grep -w "Name=$geneName" "$genomeFile" | cut -f 1)
		geneStart=$(grep -w "Name=$geneName" "$genomeFile" | cut -f 4)
		geneEnd=$(grep -w "Name=$geneName" "$genomeFile" | cut -f 5)
		#Generate image for each gene with only one gene model
		Rscript generateMergedCoverage_filtered_fromGFF3.r "$geneGff" "$outDir" "$scaffoldName" "$geneName" "$geneStart" "$geneEnd" "$inputBam1" "$inputBam2" "$inputBam3" "$inputBam4" "$inputBam5" "$inputBam6"
	else
		#Output each gene ID with multiple models to a txt file
		echo "$geneName" >> "$outIDs"
	fi

	#Clean up
	rm "$geneGff"
done < "$outDir"/tmpIDs.csv

#Clean up
rm "$outDir"/tmpIDs.csv