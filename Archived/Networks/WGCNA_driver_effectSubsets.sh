#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash WGCNA_driver_effectSubsets.sh
#Usage ex: bash WGCNA_driver_effectSubsets.sh

#Retrieve statistics outputs absolute path
#outputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")
#Retrieve analysis inputs path
#inputsPath="$outputsPath"/"$1"

#Set outputs path
#outputsPath="$outputsPath"/WGCNA_results
outputsPath="/Users/bamflappy/PfrenderLab/WGCNA_PA42_v4.1"
mkdir "$outputsPath"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputsPath directory already exsists... please remove before proceeding."
	exit 1
fi

#Retrieve factor grouping file
#grpFile="../InputData/expDesign_WGCNA_Olympics.csv"
#cp "$grpFile" "$inputsPath"

#Prepare inputs
workingDir=$outputsPath
countsTable="/Users/bamflappy/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv"
fromCol=1
toCol=24
targets="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_WGCNA_Olympics.csv"
geneCountsInter="/Users/bamflappy/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags.csv"
geneCountsTreat="/Users/bamflappy/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags.csv"
geneCountsTol="/Users/bamflappy/PfrenderLab/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags.csv"
annotIn="/Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts.csv"
annotUniprot="/Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_uniprot.csv"
annotMap="/Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_uniprotEntrezMap.csv"


#Prepare input data for WGCNA using R
#Rscript WGCNA_dataInput.r "$outputsPath" "$inputsPath"/"cleaned.csv" $2 $3 "$inputsPath"/"expDesign_Olympics.csv"
Rscript WGCNA_dataInput_effectSubsets.R $workingDir $countsTable $fromCol $toCol $targets $geneCountsInter $geneCountsTreat $geneCountsTol $annotIn $annotUniprot $annotMap $allTraits > "$outputsPath"/stats_dataInput.txt

#Construct networks
Rscript WGCNA_networkConstruction_interSubset.R $workingDir
Rscript WGCNA_networkConstruction_treatSubset.R $workingDir
Rscript WGCNA_networkConstruction_tolSubset.R $workingDir

#Relate modules to external traits
Rscript WGCNA_relateModsToExt_interSubset.R $workingDir
Rscript WGCNA_relateModsToExt_treatSubset.R $workingDir
Rscript WGCNA_relateModsToExt_tolSubset.R $workingDir

#Visualize WGCNA data
Rscript WGCNA_visualization_interSubset.R $workingDir
Rscript WGCNA_visualization_treatSubset.R $workingDir
Rscript WGCNA_visualization_tolSubset.R $workingDir

#Construct networks
Rscript WGCNA_exportNetwork_interSubset.R $workingDir
Rscript WGCNA_exportNetwork_treatSubset.R $workingDir
Rscript WGCNA_exportNetwork_tolSubset.R $workingDir

#Plot proportion of modules represented by DEG effects
#Rscript stackedBarPlot_effectSubsets_interSubset.R $workingDir
#Rscript stackedBarPlot_effectSubsets_treatSubset.R $workingDir
#Rscript stackedBarPlot_effectSubsets_tolSubset.R $workingDir
#Rscript stackedBarPlot_effectSubsets_filtered_interSubset.R $workingDir
#Rscript stackedBarPlot_effectSubsets_filtered_treatSubset.R $workingDir
#Rscript stackedBarPlot_effectSubsets_filtered_tolSubset.R $workingDir