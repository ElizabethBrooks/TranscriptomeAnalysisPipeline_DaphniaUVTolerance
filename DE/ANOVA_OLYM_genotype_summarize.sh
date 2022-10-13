#!/bin/bash

# usage: bash ANOVA_OLYM_genotype_summarize.sh

# working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"

# move to the working directory
cd $workingDir

# expression data
expData="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv"

# final summary file
sumFile="OLYM_WGCNA_aov_pairwise_summary_pValues.csv"

# add header
echo "module,treatment,genotype,interaction,Y05_Y023,Y05_E05,Y05_R2,Y023_E05,Y023_R2,E05_R2" > $sumFile

# loop over each module
counter=1
sigNum=0.05
while read line; do
	if [ $counter -gt 1 ]; then # skip header
	#if [ $counter -eq 2 ]; then # test
		# retrieve module name
		tag=$(echo "$line" | cut -d ',' -f 1)
		# set file names
		anovaSum=$(echo "${tag}_ANOVA_summary.txt" | sed "s/\"//g")
		pairSum=$(echo "${tag}_pairwise_summary.txt" | sed "s/\"//g")
		# retrieve ANOVA summary statistics
		treatment=$(cat "$anovaSum" | grep "^treatment " | tr -s '[:blank:]' ',' | cut -d ',' -f 6)				
		genotype=$(cat "$anovaSum" | grep "^genotype " | tr -s '[:blank:]' ',' | cut -d ',' -f 6)
		interaction=$(cat "$anovaSum" | grep "^genotype:treatment" | tr -s '[:blank:]' ',' | cut -d ',' -f 6)
		# retrieve pairwise summary statistics
		Y05_Y023=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "Y023" | cut -d "," -f 2)
		Y05_E05=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "E05" | cut -d "," -f 2)
		Y05_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "R2" | cut -d "," -f 2)
		Y023_E05=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,3 | grep "E05" | cut -d "," -f 2)
		Y023_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,3 | grep "R2" | cut -d "," -f 2)
		E05_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,4 | grep "R2" | cut -d "," -f 2)
		# add module results to summary file
		echo "$tag,$treatment,$genotype,$interaction,$Y05_Y023,$Y05_E05,$Y05_R2,$Y023_E05,$Y023_R2,$E05_R2" >> $sumFile
	fi
	# increment counter
	let "counter=counter+1"
done < $expData
