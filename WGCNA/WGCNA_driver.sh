#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash WGCNA_driver.sh countsFolder startColPos endColPos
#Usage ex: bash WGCNA_driver.sh sortedCoordinate_samtoolsHisat2_run3 1 24

#Retrieve statistics outputs absolute path
#outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#Retrieve analysis inputs path
#inputsPath="$outputsPath"/"$1"

#Set outputs path
#outputsPath="$outputsPath"/WGCNA_results
outputsPath="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1"
outputsPath="$outputsPath"/allGenes
mkdir "$outputsPath"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Retrieve factor grouping file
#grpFile="../InputData/expDesign_WGCNA_Olympics.csv"
#cp "$grpFile" "$inputsPath"

#Prepare input data for WGCNA using R
#Rscript WGCNA_dataInput.r "$outputsPath" "$inputsPath"/"cleaned.csv" $2 $3 "$inputsPath"/"expDesign_Olympics.csv"
Rscript WGCNA_dataInput.R

#Construct networks
Rscript WGCNA_networkConstruction.R

#Relate modules to external traits
Rscript WGCNA_relateModsToExt.R

#Visualize WGCNA data
Rscript WGCNA_visualization.R

#Construct networks
Rscript WGCNA_exportNetwork.R