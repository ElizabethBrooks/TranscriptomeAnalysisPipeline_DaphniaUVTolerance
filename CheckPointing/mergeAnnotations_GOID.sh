#Bash script to merge two GO annotation files
# and remove duplicate entries
bash mergeAnnotations_GOID.sh /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/Frozen 3.0 Genome_Annotaton/PA42_3.0_annotation_go_terms_enzyme_code_by_line_18440.txt /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/PANNZER_proteins_new/GO.out.txt

#Merge the input files
outPath=$(dirname "$2")
outFile="$outPath"/"GO_merged.txt"
cat "$1" "$2" > "$outFile"

#Remove extra tags
sed -i 's/GO://g' "$outFile"
tr -s ' ' < "$outFile" 
tr -s '\t' < "$outFile"

#Sort and remove duplicate lines
sort "$outFile" | uniq -u

#Sort by gene ID
sort -k1 -n "$outFile"

#Check number of lines
echo "Number of lines:"
wc -l "$1"
wc -l "$2"
wc -l "$outFile"

#Check number of uniqe geneIDs with GOIDs
echo "Number of gene IDs:"
genes1=$(sort -k1 -n "$1" | uniq -u | wc -l); echo "File1: $genes1"
genes2=$(sort -k1 -n "$2" | uniq -u | wc -l); echo "File2: $genes2"
genes3=$(sort -k1 -n "$outFile" | uniq -u | wc -l); echo "Merged: $genes3"