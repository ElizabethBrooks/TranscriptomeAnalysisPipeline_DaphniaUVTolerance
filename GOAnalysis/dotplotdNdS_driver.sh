#!/bin/bash

# BASH script to drive GO analysis for DE results

# usage: bash dotplotdNdS_driver.sh
# default usage ex: bash dotplotdNdS_driver.sh

# set directory for output files
inPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests"

# set outputs directory
outDir=$inPath"/GOAnalysis"

# create outputs directory
mkdir $outDir

# create a dot plot of significant GO terms for each effect set
Rscript dotplotdNdS_onlySig.R $outDir
