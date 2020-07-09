#!/bin/bash
#Bash script to merge two GO annotation files
# and remove duplicate entries
#bash mergeAnnotations_GOID_sanspanz.sh /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/Frozen3.0Genome_Annotaton/PA42_3.0_annotation_go_terms_enzyme_code_by_line_18440.txt /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/PANNZER_proteins_new/GO.out.txt

#Merge the input files
outPath=$(dirname "$2")
outFileMerged="$outPath"/"GO_merged.txt"
outFileCleaned="$outPath"/"GO_noDuplicates.txt"
cut -f 1,2 "$1" | sed 's/GO://g' > "$outFileCleaned"
tail -n +2 "$2" | cut -f 1,3 >> "$outFileCleaned"

#Remove extra tabs
unexpand -a "$outFileCleaned" > "$outFileMerged"

#Sort and remove duplicate lines
#sort "$outFileMerged" | uniq -u > "$outFileCleaned"

#Check number of lines
echo "Number of Entries"
echo "File, Total, Unique, Duplicates"
genes1a=$(wc -l "$1" | cut -d ' ' -f 1)
genes1b=$(sort "$1" | uniq -u | wc -l)
genes1c=$(($genes1a-$genes1b))
echo "File1, $genes1a, $genes1b, $genes1c"
genes2a=$(wc -l "$2" | cut -d ' ' -f 1)
genes2b=$(sort "$2" | uniq -u | wc -l)
genes2c=$(($genes2a-$genes2b))
echo "File2, $genes2a, $genes2b, $genes2c"
genes3a=$(wc -l "$outFileMerged" | cut -d ' ' -f 1)
genes3b=$(sort "$outFileMerged" | uniq -u | wc -l)
genes3c=$(($genes3a-$genes3b))
echo "Merged, $genes3a, $genes3b, $genes3c"

#Check number of geneIDs
printf "\nNumber of Gene IDs\n"
echo "File, Total, Unique, Duplicates"
genes1a=$(cut -f 1 "$1" | wc -l)
genes1b=$(cut -f 1 "$1" | uniq -u | wc -l)
genes1c=$(($genes1a-$genes1b))
echo "File1, $genes1a, $genes1b, $genes1c"
genes2a=$(cut -f 1 "$2" | wc -l)
genes2b=$(cut -f 1 "$2" | uniq -u | wc -l)
genes2c=$(($genes2a-$genes2b))
echo "File2, $genes2a, $genes2b, $genes2c"
genes3a=$(cut -f 1 "$outFileMerged" | wc -l)
genes3b=$(cut -f 1 "$outFileCleaned" | uniq -u | wc -l)
genes3c=$(($genes3a-$genes3b))
echo "Merged, $genes3a, $genes3b, $genes3c"

#Check number of GOIDs
printf "\nNumber of GO IDs\n"
echo "File, Total, Unique, Duplicates"
genes1a=$(cut -f 2 "$1" | wc -l)
genes1b=$(cut -f 2 "$1" | uniq -u | wc -l)
genes3c=$(($genes1a-$genes1b))
echo "File1, $genes1a, $genes1b, $genes1c"
genes2a=$(cut -f 2 "$2" | wc -l)
genes2b=$(cut -f 2 "$2" | uniq -u | wc -l)
genes3c=$(($genes2a-$genes2b))
echo "File2, $genes2a, $genes2b, $genes2c"
genes3a=$(cut -f 2 "$outFileMerged" | wc -l)
genes3b=$(cut -f 2 "$outFileCleaned" | uniq -u | wc -l)
genes3c=$(($genes3a-$genes3b))
echo "Merged, $genes3a, $genes3b, $genes3c"
