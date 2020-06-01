#Bash script to merge two GO annotation files
# and remove duplicate entries
#bash mergeAnnotations_GOID.sh /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/Frozen3.0Genome_Annotaton/PA42_3.0_annotation_go_terms_enzyme_code_by_line_18440.txt /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/PANNZER_proteins_new/GO.out.txt

#Merge the input files
outPath=$(dirname "$2")
outFileMerged="$outPath"/"GO_merged.txt"
outFileCleaned="$outPath"/"GO_noDuplicates.txt"
cut -f 1,2 "$1" | sed 's/GO://g' > "$outFileCleaned"
tail -n +2 "$2" | cut -f 1,3 >> "$outFileCleaned"

#Remove extra tabs
unexpand -a "$outFileCleaned" > "$outFileMerged"

#Sort and remove duplicate lines
sort "$outFileMerged" | uniq -u > "$outFileCleaned"

#Check number of lines
echo "Number of lines..."
genes1=$(wc -l "$1" | cut -d ' ' -f 1); echo "File1: $genes1"
genes2=$(wc -l "$2" | cut -d ' ' -f 1); echo "File2: $genes2"
genes3=$(wc -l "$outFileMerged" | cut -d ' ' -f 1); echo "Merged: $genes3"
genes4=$(wc -l "$outFileCleaned" | cut -d ' ' -f 1); echo "No Duplicates: $genes4"
genes5=$(($genes3-$genes4)); echo "Number of Duplicates: $genes5"

#Check number of uniqe geneIDs with GOIDs
printf "\nNumber of gene IDs...\n"
genes1=$(sort -k1 -n "$1" | uniq -u | wc -l); echo "File1: $genes1"
genes2=$(sort -k1 -n "$2" | uniq -u | wc -l); echo "File2: $genes2"
genes3=$(sort -k1 -n "$outFileMerged" | uniq -u | wc -l); echo "Merged: $genes3"
genes4=$(sort -k1 -n "$outFileCleaned" | uniq -u | wc -l); echo "No Duplicates: $genes4"
genes5=$(($genes3-$genes4)); echo "Number of Duplicates: $genes5"