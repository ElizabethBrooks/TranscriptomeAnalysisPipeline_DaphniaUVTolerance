#!/bin/bash

# usage: bash ANOVA_OLYM_tolerance_aov_driver.sh

# set tag
setTag="OLYM_60"

# working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA"
workingDir=$workingDir"/"$setTag

# create the working directory
mkdir $workingDir

# design file
designFile="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsTolerance.csv"
#designFile="/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/expDesign_OlympicsGenotypes.csv"

# expression data
expData="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/ensembl/GCA_021134715.1/biostatistics/NetworkAnalysis/WGCN_tolerance_WGCNA/OLYM_60_eigengeneExpression.csv"

# loop over each module
counter=1
while read line; do
	if [ $counter -eq 1 ]; then
		echo "$line" > $workingDir/"tmpHeader.csv"
	else
		# setup expression data
		cat $workingDir/"tmpHeader.csv" > $workingDir/"tmpMod.csv"
		echo "$line" >> $workingDir/"tmpMod.csv"
		# calcluate ANOVA
		Rscript ANOVA_OLYM_tolerance_aov.r $workingDir $workingDir/"tmpMod.csv" $designFile
		#Rscript ANOVA_OLYM_genotype_aov.r $workingDir $workingDir/"tmpMod.csv" $designFile
		# clean up
		rm $workingDir/"tmpMod.csv"
	fi
	let "counter=counter+1"
done < $expData

# clean up
rm $workingDir/"tmpHeader.csv"