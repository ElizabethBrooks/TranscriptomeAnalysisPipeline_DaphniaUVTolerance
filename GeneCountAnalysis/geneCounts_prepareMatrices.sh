#!/bin/bash
#Subset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Legacy gene counts were created using trimmomatic, tophat, and cuffdiff

#Directory for gene count tables
prefix="../../GeneCounts_Merged/"

#Merge counts from subset of new gene count table matching samples
# in the previously created "legacy" gene count table
sed -i -e 's|^|stats_hisat2EdgeR_test1/stats_hisat2EdgeR_counts_test1/|' "$prefix"mergeGuideFile.txt
python "$prefix"merge_tables.py "$prefix"mergeGuideFile.txt

#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$prefix"merged_counts_subset.txt | tr -s '[:blank:]' ',' > "$prefix"merged_counts_subset_cleaned.csv
cat "$prefix"merged_counts_legacy.txt | tr -s '[:blank:]' ',' > "$prefix"merged_counts_legacy_cleaned.csv

#Remove extra lines from the new gene count table
egrep -v "unique|ambiguous|feature|aligned|aQual" "$prefix"merged_counts_subset_cleaned.csv > "$prefix"merged_counts_subset_rowCleaned.csv

#Retrieve the sorted row tags of gene IDs for each file
cat "$prefix"merged_counts_subset_rowCleaned.csv | cut -d, -f1 | sort -d > "$prefix"merged_counts_subset_rowCleanedTags.csv
cat "$prefix"merged_counts_legacy_cleaned.csv | cut -d, -f1 | sort -d > "$prefix"merged_counts_legacy_rowTags.csv

#Compare row tags to identify lines unique to either file,
# in order to ensure matching gene IDs before merging
#Supress output of row tags found in both files using the -3 flag
comm -3 "$prefix"merged_counts_subset_rowCleanedTags.csv "$prefix"merged_counts_legacy_rowTags.csv > "$prefix"merged_counts_uniqueRowTags.csv

#Add tags to each sample name in each file indicating
# the alignment method used, which is the primary difference
# between the gene count tables beinf compared (subset and legacy)
sed 's/Pool1/Pool1_H/g' "$prefix"merged_counts_subset_rowCleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$prefix"merged_counts_subset_tagged.csv
sed 's/Pool1/Pool1_T/g' "$prefix"merged_counts_legacy_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$prefix"merged_counts_legacy_tagged.csv

#Remove the row tags with gene IDs before merging
sed 's/\([^,]*\),\(.*\)/\2/' "$prefix"merged_counts_subset_tagged.csv > "$prefix"merged_counts_subset_trimmed.csv
#sed 's/\([^,]*\),\(.*\)/\2/' "$prefix"merged_counts_legacy_tagged.csv > "$prefix"merged_counts_legacy_trimmed.csv

#Merge both the subset and legacy gene count tables for further comparison
paste -d , "$prefix"merged_counts_legacy_tagged.csv merged_counts_subset_noIDs.csv > "$prefix"final_merged_counts.csv

#Trim off the row tags with gene IDs from the final merged file
sed 's/\([^,]*\),\(.*\)/\2/' "$prefix"final_merged_counts.csv > "$prefix"final_merged_counts_trimmed.csv