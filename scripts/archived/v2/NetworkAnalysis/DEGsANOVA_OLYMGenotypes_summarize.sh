#!/bin/bash

# usage: bash DEGsANOVA_OLYMGenotypes_summarize.sh analysisType set
# usage ex: bash DEGsANOVA_OLYMGenotypes_summarize.sh Tolerance OLYM

# working directory
workingDir=$1

# retrieve analysis type
analysisType=$2

# retrieve set
set=$3

# minimum module size
minModSize=$4

#Create directory for output files
inDir=$(grep "WGCNA:" ../InputData/outputPaths.txt | tr -d " " | sed "s/WGCNA://g")
inDir=$inDir"/"$analysisType

# move to the working directory
cd $workingDir

# expression data
expData=$inDir"/"$set"_"$minModSize"_eigengeneExpression.csv"

# final summary files
sumFileAov="aov_summary_pValues.csv"

# add headers
echo "module,treatment,tolerance,toleranceGenotype,treatmentTolerance,treatmentToleranceGenotype" > $sumFileAov

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
		treatment=$(cat "$anovaSum" | grep "^treatment " | tr -s '[:blank:]' ',' | sed "s/<//g" | sed "s/,,/,/g" | cut -d ',' -f 6 | head -1)	
		tolerance=$(cat "$anovaSum" | grep "^tolerance " | tr -s '[:blank:]' ',' | sed "s/<//g" | sed "s/,,/,/g" | cut -d ',' -f 6 | head -1)
		toleranceGenotype=$(cat "$anovaSum" | grep "^tolerance:genotype " | tr -s '[:blank:]' ',' | sed "s/<//g" | sed "s/,,/,/g" | cut -d ',' -f 6 | head -1)			
		treatmentTolerance=$(cat "$anovaSum" | grep "^treatment:tolerance " | tr -s '[:blank:]' ',' | sed "s/<//g" | sed "s/,,/,/g" | cut -d ',' -f 6 | head -1)
		treatmentToleranceGenotype=$(cat "$anovaSum" | grep "^treatment:tolerance:genotype " | tr -s '[:blank:]' ',' | sed "s/<//g" | sed "s/,,/,/g" | cut -d ',' -f 6 | head -1)
		# add module results to summary files
		echo "$tag,$treatment,$tolerance,$toleranceGenotype,$treatmentTolerance,$treatmentToleranceGenotype" >> $sumFileAov
	fi
	# increment counter
	let "counter=counter+1"
done < $expData
