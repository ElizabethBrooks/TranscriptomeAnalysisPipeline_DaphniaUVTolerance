#!/bin/bash

path=$(echo /Users/bamflappy/PfrenderLab/PA42_v4.1/genomicData_PA42_v4.1)

grep ">" "$path"/PA42.4.1.pep_longest.fasta | sed "s/-pep//g" | sed "s/>//g" > "$path"/col1.txt
tr -d "\n\r" < "$path"/col1.txt | sed "s/dp\_gene/\ndp\_gene/g" > "$path"/col1_cleaned.txt

grep ">" "$path"/PA42.4.1.pep_longest.fasta | sed "s/>//g" > "$path"/col2.txt
tr -d "\n\r" < "$path"/col2.txt | sed "s/dp\_gene/\ndp\_gene/g" > "$path"/col2_cleaned.txt

paste "$path"/col1_cleaned.txt "$path"/col2_cleaned.txt > "$path"/cols_merged.txt

tail -n +2 "$path"/cols_merged.txt > "$path"/PA42.4.1.fasta.gene_pep_map

#Clean up 
rm "$path"/col*.txt