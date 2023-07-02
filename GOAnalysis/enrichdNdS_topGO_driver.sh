#!/bin/bash

# BASH script to drive GO analysis for dN/dS results

# usage: bash enrichdNdS_topGO_driver.sh
# default usage ex: bash enrichdNdS_topGO_driver.sh

# retrieve analysis type
analysisType=$1

# set directory for output files
inPath="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/selectionTests"

# retrieve dN dS values
inFile=$inPath"/Pulex_Olympics_kaksResults_dNdS_cleaned.csv"

# retrieve functional annotations
GOmaps=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")
GOmaps=$(dirname $GOmaps)
GOmaps=$GOmaps"/geneToGO_tagged_map.txt"

# set outputs directory
outDir=$inPath"/GOAnalysis"

# create outputs directory
mkdir $outDir

# determine the direction of expression for genes under negative selection
Rscript enrichdNdS_negative_topGO.R $inFile $GOmaps $outDir

# determine the direction of expression for genes under positive selection
Rscript enrichdNdS_positive_topGO.R $inFile $GOmaps $outDir

# determine the direction of expression for genes under positive selection
Rscript enrichdNdS_strict_topGO.R $inFile $GOmaps $outDir

# determine the direction of expression for genes under positive selection
Rscript enrichdNdS_error_topGO.R $inFile $GOmaps $outDir
