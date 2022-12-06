#!/bin/bash
#Script to run Rscripts that perform WGCNA analysis of OLYM data
#Usage: bash WGCNA_driver_genotype.sh
#Usage ex: bash WGCNA_driver_genotype.sh

# input counts file
inputCounts="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts_logTransformed.csv"

# experimental design file
expDesign="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv"

# name output directory
outDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics"

# create output directory
mkdir $outDir

# WGCNA data input
Rscript dataInput_consensus_genotype_WGCNA.R $outDir $inputCounts 1 24 $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 1 24 OLYM $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 1 12 Tol $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 13 24 NTol $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 1 6 Y05 $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 7 12 Y023 $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 13 18 E05 $expDesign
Rscript dataInput_set_WGCNA.R $outDir $inputCounts 1 24 19 24 R2 $expDesign

# WGCNA soft power picking
Rscript pickSoftPower_consensus_WGCNA.R $outDir genotype
Rscript pickSoftPower_set_WGCNA.R $outDir OLYM
Rscript pickSoftPower_set_WGCNA.R $outDir Tol
Rscript pickSoftPower_set_WGCNA.R $outDir NTol
Rscript pickSoftPower_set_WGCNA.R $outDir Y05
Rscript pickSoftPower_set_WGCNA.R $outDir Y023
Rscript pickSoftPower_set_WGCNA.R $outDir E05
Rscript pickSoftPower_set_WGCNA.R $outDir R2

# WGCNA network construction
Rscript networkConstruction_consensus_genotype_WGCNA.R $outDir 14 60
Rscript networkConstruction_set_WGCNA.R $outDir OLYM 8 60
Rscript networkConstruction_set_WGCNA.R $outDir Tol 14 60
Rscript networkConstruction_set_WGCNA.R $outDir NTol 14 60
Rscript networkConstruction_set_WGCNA.R $outDir Y05 20 60
Rscript networkConstruction_set_WGCNA.R $outDir Y023 12 60
Rscript networkConstruction_set_WGCNA.R $outDir E05 14 60
Rscript networkConstruction_set_WGCNA.R $outDir R2 9 60

# WGCNA consensus network analsysis
Rscript networkAnalysis_consensus_WGCNA.R $outDir OLYM 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir Tol 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir NTol 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir Y05 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir Y023 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir E05 60 genotype
Rscript networkAnalysis_consensus_WGCNA.R $outDir R2 60 genotype

# WGCNA compare consensus eigengene networks
Rscript eigengeneNetworks_consensus_genotype_WGCNA.R $outDir 60

# WGCNA export eigengene expression values
Rscript eigengeneExpression_set_WGCNA.R $outDir OLYM 60
Rscript eigengeneExpression_set_WGCNA.R $outDir Tol 60
Rscript eigengeneExpression_set_WGCNA.R $outDir NTol 60
Rscript eigengeneExpression_set_WGCNA.R $outDir Y05 60
Rscript eigengeneExpression_set_WGCNA.R $outDir Y023 60
Rscript eigengeneExpression_set_WGCNA.R $outDir E05 60
Rscript eigengeneExpression_set_WGCNA.R $outDir R2 60
