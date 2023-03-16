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

# set filter type
filterType="filteredMapQ"

# run script to add read groups
bash addReadGroups_samtools.sh $sortedFolderName

# run script to filter bam files by mapq
bash filterByMapQ_samtools.sh $sortedFolderName $filterType

# run variant calling script
bash variantCallingMerged_bcftools.sh $sortedFolderName $filterType

# run variant filtering script
bash variantFiltering_bcftools.sh $sortedFolderName $filterType

# run script to generate consensus sequences
bash generateConsensusMerged_bcftools.sh $sortedFolderName $filterType
