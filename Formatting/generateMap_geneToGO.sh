#!/bin/bash

# script to map gene IDs with GO annotations
# usage: bash generateMap_geneToGO.sh

# retrieve functional annotations
annotations=$(grep "functionalAnnotations:" ../InputData/inputPaths.txt | tr -d " " | sed "s/functionalAnnotations://g")

# name output directory for the analysis type
outDir=$(dirname $annotations)

# move to outputs directory
cd $outDir

# retrieve gene and GO annotations
cat $annotations | cut -f 1,9 | grep "ontology_id" | sed 's/ontology_id\ //g' | sed 's/\"//g' > "geneToGO_map.txt"

# add gene- tag
cat "geneToGO_map.txt" | sed 's/^/gene-/g' | sed 's/\"//g' > "geneToGO_tagged_map.txt"
