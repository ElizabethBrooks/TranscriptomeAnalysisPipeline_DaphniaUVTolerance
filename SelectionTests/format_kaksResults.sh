#!/bin/bash

# script to format results from tests for selection
# usage: bash format_kaksResults.sh

# retrieve current working directory
currDir=$(pwd)

# retrieve base directory path
baseDir=$(dirname $currDir)

# set software path
softwarePath=$(grep "pal2nal:" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/pal2nal://g")

# retrieve genome features absolute path for alignment
genomeFeatures=$(grep "genomeFeatures" $baseDir"/InputData/softwarePaths.txt" | tr -d " " | sed "s/genomeFeatures://g")

# retrieve sorted reads input absolute path
inputsPath=$(grep "aligningGenome:" $baseDir"/InputData/outputPaths.txt" | tr -d " " | sed "s/aligningGenome://g")

# retrieve genome reference absolute path for alignment
refPath=$(grep "genomeReference" $baseDir"/InputData/inputPaths.txt" | tr -d " " | sed "s/genomeReference://g")

# set inputs path
inputsPath=$inputsPath"/variantsCalled_samtoolsBcftools"

# retrieve input bam file type
type="filteredMapQ"

# make outputs directory name
outFolder=$inputsPath"/selectionTests"

# set inputs folder
inputsPath=$inputsPath"/features_gffread"

# retrieve file name of reference
refTag=$(basename $refPath)

# set results file path
resultsFile=$outFolder"/kaksResults.csv"

# status message
echo "Formatting results..."

# merge each of the ka ks results files
tail -n+2 $outFolder"/kaksResults_"*".csv" >> $resultsFile

# fix formatting of the results file
finalResults="$outFolder"/Pulex_Olympics_kaksResults.csv
cat $resultsFile | grep -v "^==>" | tr -s '[:blank:]' ',' | sed "s/=,/=/g" | sed "s/dN\/dS=//g" | sed "s/dN,=//g" | sed "s/dS,=//g" | sed "s/t=//g" | sed "s/S=//g" | sed "s/N=//g" > $finalResults

# remove genes with no dNdS value
finalResults="$outFolder"/Pulex_Olympics_kaksResults_dNdS_cleaned.csv
cat $resultsFile | grep -v "^==>" | awk '$4!=""' | tr -s '[:blank:]' ',' | sed "s/=,/=/g" | sed "s/dN\/dS=//g" | sed "s/dN,=//g" | sed "s/dS,=//g" | sed "s/t=//g" | sed "s/S=//g" | sed "s/N=//g" > $finalResults

# keep only genes with dN/dS= 99.0000
finalResults="$outFolder"/Pulex_Olympics_kaksResults_dNdS_99.csv
echo "geneID,t,S,N,dNdS,dN,dS" > $finalResults
cat $resultsFile | grep "dN/dS= 99.0000" | tr -s '[:blank:]' ',' | sed "s/=,/=/g" | sed "s/dN\/dS=//g" | sed "s/dN,=//g" | sed "s/dS,=//g" | sed "s/t=//g" | sed "s/S=//g" | sed "s/N=//g" >> $finalResults

# keep only genes without dNdS values
finalResults="$outFolder"/Pulex_Olympics_kaksResults_dNdS_NA.csv
cat $resultsFile | awk '$4==""' | tr -s '[:blank:]' ',' | sed "s/=,/=/g" | sed "s/dN\/dS=//g" | sed "s/dN,=//g" | sed "s/dS,=//g" | sed "s/t=//g" | sed "s/S=//g" | sed "s/N=//g" > $finalResults

# clean up
rm $outFolder"/kaksResults_"*".csv"
rm $resultsFile

# status message
echo "Analysis complete!"
