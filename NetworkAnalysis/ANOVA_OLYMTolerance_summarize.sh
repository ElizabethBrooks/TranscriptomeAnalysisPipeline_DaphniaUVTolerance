#!/bin/bash

# usage: bash ANOVA_OLYMTolerance_summarize.sh analysisType set
# usage ex: bash ANOVA_OLYMTolerance_summarize.sh Tolerance OLYM

# retrieve analysis type
analysisType=$1

#Create directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# working directory
workingDir=$inDir"/ANOVA"

# move to the working directory
cd $workingDir

# retrieve set
set=$2

# minimum module size
minModSize=30

# expression data
expData=$inDir"/"$set"_"$minModSize"_eigengeneExpression.csv"

# final summary files
sumFileAov="aov_summary_pValues.csv"

# add headers
echo "module,treatment,tolerance,interaction" > $sumFileAov

# loop over each module
counter=1
sigNum=0.05
while read line; do
	if [ $counter -gt 1 ]; then # skip header
	#if [ $counter -eq 2 ]; then # test
		# retrieve module name
		tag=$(echo "$line" | cut -d ',' -f 1)
		# set file names
		anovaSum=$(echo "${tag}_ANOVA_summary.csv" | sed "s/\"//g")
		tukeySum=$(echo "${tag}_tukey_summary.txt" | sed "s/\"//g")
		pairSum=$(echo "${tag}_pairwise_summary.txt" | sed "s/\"//g")
		# retrieve ANOVA summary statistics
		treatment=$(cat "$anovaSum" | grep "^treatment " | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")				
		tolerance=$(cat "$anovaSum" | grep "^tolerance " | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")
		interaction=$(cat "$anovaSum" | grep "^tolerance:treatment" | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")
		# add module results to summary files
		echo "$tag,$treatment,$tolerance,$interaction" >> $sumFileAov
	fi
	# increment counter
	let "counter=counter+1"
done < $expData
