#!/bin/bash

# usage: bash ANOVA_OLYM_genotype_summarize.sh

# working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCNA_DEGenotypes"

# move to the working directory
cd $workingDir

# expression data
expData="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/WGCN_OLYM_WGCNA/OLYM_60_eigengeneExpression.csv"

# final summary files
sumFileAov="OLYM_WGCNA_aov_tukey_summary_pValues.csv"
sumFilePair="OLYM_WGCNA_pairwise_summary_pValues.csv"

# add headers
echo "module,treatment,genotype,interaction,Y023_Y05,E05_Y05,R2_Y05,E05_Y023,R2_Y023,R2_E05" > $sumFileAov
echo "module,Y05_Y023,Y05_E05,Y05_R2,Y023_E05,Y023_R2,E05_R2" > $sumFilePair

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
		tukeySum=$(echo "${tag}_tukey_summary.txt" | sed "s/\"//g")
		pairSum=$(echo "${tag}_pairwise_summary.txt" | sed "s/\"//g")
		# retrieve ANOVA summary statistics
		treatment=$(cat "$anovaSum" | grep "^treatment " | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")				
		genotype=$(cat "$anovaSum" | grep "^genotype " | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")
		interaction=$(cat "$anovaSum" | grep "^genotype:treatment" | tr -s '[:blank:]' ',' | cut -d ',' -f 6 | sed "s/<//g")
		# retrieve tukey post hoc test statistics
		Y023_Y05=$(cat $tukeySum | grep "Y023-Y05" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		E05_Y05=$(cat $tukeySum | grep "E05-Y05" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		R2_Y05=$(cat $tukeySum | grep "R2-Y05" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		E05_Y023=$(cat $tukeySum | grep "E05-Y023" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		R2_Y023=$(cat $tukeySum | grep "R2-Y023" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		R2_E05=$(cat $tukeySum | grep "R2-E05" | tr -s '[:blank:]' ',' | cut -d "," -f 5 | sed "s/<//g")
		# retrieve pairwise summary statistics
		Y05_Y023=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "Y023" | cut -d "," -f 2 | sed "s/<//g")
		Y05_E05=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "E05" | cut -d "," -f 2 | sed "s/<//g")
		Y05_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,2 | grep "R2" | cut -d "," -f 2 | sed "s/<//g")
		Y023_E05=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,3 | grep "E05" | cut -d "," -f 2 | sed "s/<//g")
		Y023_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,3 | grep "R2" | cut -d "," -f 2 | sed "s/<//g")
		E05_R2=$(cat $pairSum | tail -5 | head -3 | tr -s '[:blank:]' ',' | cut -d "," -f 1,4 | grep "R2" | cut -d "," -f 2 | sed "s/<//g")
		# add module results to summary files
		echo "$tag,$treatment,$genotype,$interaction,$Y023_Y05,$E05_Y05,$R2_Y05,$E05_Y023,$R2_Y023,$R2_E05" >> $sumFileAov
		echo "$tag,$Y05_Y023,$Y05_E05,$Y05_R2,$Y023_E05,$Y023_R2,$E05_R2" >> $sumFilePair
	fi
	# increment counter
	let "counter=counter+1"
done < $expData
