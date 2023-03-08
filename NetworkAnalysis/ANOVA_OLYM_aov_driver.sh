#!/bin/bash

# usage: bash ANOVA_OLYM_aov_driver.sh analysisType set
# default usage ex: bash ANOVA_OLYM_aov_driver.sh Tolerance OLYM
# usage ex: bash ANOVA_OLYM_aov_driver.sh Tolerance Tol
# usage ex: bash ANOVA_OLYM_aov_driver.sh Tolerance NTol

# retrieve analysis type
analysisType=$1

# retrieve set
set=$2

# minimum module size
minModSize=30

#Create directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# working directory
workingDir=$inDir"/ANOVA_"$set"_"$minModSize

# create the working directory
mkdir $workingDir

# get current directory
currDir=$(pwd)

# retrieve experimental design data path
cd $(dirname "../InputData/expDesign_Olympics"$analysisType".csv")
designPath=$(pwd)
designFile=$designPath"/expDesign_Olympics"$analysisType".csv"

# move back to scripts directory
cd $currDir

# expression data
expData=$inDir"/"$set"_"$minModSize"_eigengeneExpression.csv"

# loop over each module
counter=1
while read line; do
	if [ $counter -eq 1 ]; then
		echo "$line" > $workingDir"/tmpHeader.csv"
	else
		# setup expression data
		cat $workingDir"/tmpHeader.csv" > $workingDir"/tmpMod.csv"
		echo "$line" >> $workingDir"/tmpMod.csv"
		# calcluate ANOVA
		Rscript "ANOVA_OLYM"$analysisType"_aov.R" $workingDir $workingDir"/tmpMod.csv" $designFile
		# clean up
		rm $workingDir"/tmpMod.csv"
	fi
	let "counter=counter+1"
done < $expData

# clean up
rm $workingDir"/tmpHeader.csv"

# summarize the ANOVA results
bash ANOVA_OLYMTolerance_summarize.sh $analysisType $set