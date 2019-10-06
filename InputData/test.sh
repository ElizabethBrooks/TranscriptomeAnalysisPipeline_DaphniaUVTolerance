#!/bin/bash
#Retrieve inputs for gff absolute path
inputsFile="statsInputs_edgeR.txt"
genomeFile=$(head -n 1 $inputsFile)
echo "GENOME FILE: $genomeFile"