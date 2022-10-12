#!/bin/bash
#Script to run Rscripts that perform WGCNA analysis of OLYM data
#Usage: bash WGCNA_driver_OLYM_trimmed.sh
#Usage ex: bash WGCNA_driver_OLYM_trimmed.sh

# WGCNA data input
Rscript dataInput_set_trimmed_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv 1 24 1 24 OLYM /Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_tolerance_WGCNA_Olympics.csv

# WGCNA soft power picking
Rscript pickSoftPower_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM

## minModuleSize = 30
# WGCNA network construction
Rscript networkConstruction_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 8 30
# WGCNA export eigengene expression values
Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 30

## minModuleSize = 60
# WGCNA network construction
Rscript networkConstruction_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 8 60
# WGCNA export eigengene expression values
Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 60

## minModuleSize = 90
# WGCNA network construction
Rscript networkConstruction_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 8 90
# WGCNA export eigengene expression values
Rscript eigengeneExpression_set_WGCNA.R /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_trimmed_WGCNA OLYM 90