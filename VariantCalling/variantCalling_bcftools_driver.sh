#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -pe smp 4
#$ -N VC_pipeline_jobOutput

# script to perform bam read quaity filtering
# usage: qsub variantCalling_bcftools_driver.sh sortedFolderName
# usage Ex: qsub variantCalling_bcftools_driver.sh sortedCoordinate_samtoolsHisat2_run1

#Required modules for ND CRC servers
module load bio

# set sorted folder name
sortedFolderName="$1"

# retrieve aligned reads input absolute path
inputsPath=$(grep "aligningGenome:" ../InputData/outputPaths.txt | tr -d " " | sed "s/aligningGenome://g")

# set sorting outputs absolute path
inputsDir=$inputsPath"/"$sortedFolderName

# set filter type
filterType="filteredMapQ"

# merge alignments for each genotype
bash sortingMerge_samtools.sh $sortedFolderName E05
bash sortingMerge_samtools.sh $sortedFolderName R2
bash sortingMerge_samtools.sh $sortedFolderName Y05
bash sortingMerge_samtools.sh $sortedFolderName Y023

# run script to add read groups
bash addReadGroups_samtools.sh

# run script to filter bam files by mapq
bash filterByMapQ_samtools.sh

# run variant calling script
bash variantCallingMerged_bcftools.sh

# run variant filtering script
bash variantFiltering_bcftools.sh

# run script to generate consensus sequences
bash generateConsensusMerged_bcftools.sh
bash generateConsensusGenotype_bcftools.sh "E05"
bash generateConsensusGenotype_bcftools.sh "R2"
bash generateConsensusGenotype_bcftools.sh "Y05"
bash generateConsensusGenotype_bcftools.sh "Y023"

#Copy previous summaries
cp "$inputsDir"/*.txt "$outputFolder"

# clean up
#rm -r $inputsDir
