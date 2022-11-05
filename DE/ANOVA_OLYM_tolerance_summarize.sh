#!/bin/bash

# usage: bash ANOVA_OLYM_tolerance_summarize.sh

# set tag
setTag="OLYM_60"

# working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA"
workingDir=$workingDir"/"$setTag

# move to the working directory
cd $workingDir

# expression data
expData="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA/OLYM_60_eigengeneExpression.csv"

# final summary files
sumFileAov="OLYM_WGCNA_aov_summary_pValues.csv"

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
