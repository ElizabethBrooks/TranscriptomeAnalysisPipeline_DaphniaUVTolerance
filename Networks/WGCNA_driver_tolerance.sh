#!/bin/bash
#Script to run Rscripts that perform WGCNA analysis of OLYM data
#Usage: bash WGCNA_driver_tolerance.sh
#Usage ex: bash WGCNA_driver_tolerance.sh

# WGCNA data input
Rscript dataInput_consensus_tolerance_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 12 tol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 24 nTol /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 6 Y05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 7 12 Y023 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 13 18 E05 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv
Rscript dataInput_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 19 24 R2 /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_treatment_WGCNA_Olympics.csv

# WGCNA soft power picking
Rscript pickSoftPower_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA tolerance
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA tol
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA nTol
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y05
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y023
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA E05
Rscript pickSoftPower_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA R2

# WGCNA network construction
Rscript networkConstruction_consensus_tolerance_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA 14 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA tol 14 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA nTol 14 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y05 20 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y023 12 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA E05 14 60
Rscript networkConstruction_subset_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA R2 9 60

# WGCNA consensus network analsysis
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA tol 60 tolerance
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA nTol 60 tolerance
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y05 60 tolerance
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA Y023 60 tolerance
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA E05 60 tolerance
Rscript networkAnalysis_consensus_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_tolerance_WGCNA R2 60 tolerance
