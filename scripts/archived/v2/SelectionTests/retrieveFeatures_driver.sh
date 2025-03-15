#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N retrieveFeatures_jobOutput

# script to retrieve features for tests for selection
# usage: qsub retrieveFeatures_driver.sh

#Required modules for ND CRC servers
module load bio/2.0

# retrieve features for the reference and consensus genomes
bash retrieveFeaturesMerged_gffread.sh
#bash retrieveFeaturesGenotype_gffread.sh

# retreive protein coding sequence transcript names
bash retrieve_proteinCoding.sh
