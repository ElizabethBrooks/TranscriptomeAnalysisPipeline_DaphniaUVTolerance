#!/bin/bash

# usage: bash ANOVA_OLYM_genotype_summarize.sh

# working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"

# expression data
expData="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv"

# loop over each module
counter=1
while read line; do
	if [ $counter -gt 1 ]; then # skip header
		# retrieve module name
		tag=$(echo "$line" | cut -d ',' -f 1)
		echo $tag
		# retrieve summary statistics

	fi
	let "counter=counter+1"
done < $expData
