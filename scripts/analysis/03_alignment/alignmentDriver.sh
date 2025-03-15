#!/bin/bash

# script to perform alignment of sequencing data to a reference
# usage: bash alignmentDriver.sh trimmedFolder
# usage Ex: bash alignmentDriver.sh trimmed

# retrieve input folder of trimmed data
inputFolder=$1

# run script to build the reference
bash building_hisat2.sh $inputFolder

# run slected alignment software
bash alignment_hisat2.sh $inputFolder
