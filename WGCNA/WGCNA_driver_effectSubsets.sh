#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash WGCNA_driver_effectSubsets.sh countsFolder
#Usage ex: bash WGCNA_driver_effectSubsets.sh sortedCoordinate_samtoolsHisat2_run3

#Retrieve statistics outputs absolute path
#outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#Retrieve analysis inputs path
#inputsPath="$outputsPath"/"$1"

#Set outputs path
#outputsPath="$outputsPath"/WGCNA_results
outputsPath="/home/mae/Documents/RNASeq_Workshop_ND/WGCNA_PA42_v4.1"
outputsPath="$outputsPath"/effectSubsets
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
Rscript WGCNA_dataInput_effectSubsets.R > "$outputsPath"/stats_dataInput.txt

#Construct networks
Rscript WGCNA_networkConstruction_interSubset.R
Rscript WGCNA_networkConstruction_treatSubset.R
Rscript WGCNA_networkConstruction_tolSubset.R

#Relate modules to external traits
Rscript WGCNA_relateModsToExt_interSubset.R
Rscript WGCNA_relateModsToExt_treatSubset.R
Rscript WGCNA_relateModsToExt_tolSubset.R

#Visualize WGCNA data
Rscript WGCNA_visualization_interSubset.R
Rscript WGCNA_visualization_treatSubset.R
Rscript WGCNA_visualization_tolSubset.R

#Construct networks
Rscript WGCNA_exportNetwork_interSubset.R
Rscript WGCNA_exportNetwork_treatSubset.R
Rscript WGCNA_exportNetwork_tolSubset.R

#Plot proportion of modules represented by DEG effects
Rscript stackedBarPlot_effectSubsets_interSubset.R
Rscript stackedBarPlot_effectSubsets_treatSubset.R
Rscript stackedBarPlot_effectSubsets_tolSubset.R
Rscript stackedBarPlot_effectSubsets_filtered_interSubset.R
Rscript stackedBarPlot_effectSubsets_filtered_treatSubset.R
Rscript stackedBarPlot_effectSubsets_filtered_tolSubset.R