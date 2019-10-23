#!/bin/bash
#Subset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Legacy gene counts were created using trimmomatic, tophat, and cuffdiff

#Merge counts from subset of new gene count table matching samples
# in the previously created "legacy" gene count table
sed -i -e 's|^|stats_hisat2EdgeR_test1/stats_hisat2EdgeR_counts_test1/|' mergeGuideFile.txt
python merge_tables.py mergeGuideFile.txt

#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat merged_counts_subset.txt | tr -s '[:blank:]' ',' > merged_counts_subset_cleaned.csv
cat merged_counts_legacy.txt | tr -s '[:blank:]' ',' > merged_counts_legacy_cleaned.csv

#Remove extra lines from the new gene count table
egrep -v "unique|ambiguous|feature|aligned|aQual" merged_counts_subset_cleaned.csv > merged_counts_subset_rowCleaned.csv

#Retrieve the sorted row tags of gene IDs for each file
cat merged_counts_subset_rowCleaned.csv | cut -d, -f1 | sort -d > merged_counts_subset_rowCleanedTags.csv
cat merged_counts_legacy_cleaned.csv | cut -d, -f1 | sort -d > merged_counts_legacy_rowTags.csv

#Compare row tags to identify lines unique to either file,
# in order to ensure matching gene IDs before merging
#Supress output of row tags found in both files using the -3 flag
comm -3 merged_counts_subset_rowCleanedTags.csv merged_counts_legacy_rowTags.csv > merged_counts_uniqueRowTags.csv

#Add tags to each sample name in each file indicating
# the alignment method used, which is the primary difference
# between the gene count tables beinf compared (subset and legacy)
sed 's/Pool1/Pool1_H/g' merged_counts_subset_rowCleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > merged_counts_subset_tagged.csv
sed 's/Pool1/Pool1_T/g' merged_counts_legacy_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > merged_counts_legacy_tagged.csv

#Remove the row tags with gene IDs from the subset file before merging
sed 's/\([^,]*\),\(.*\)/\2/' merged_counts_subset_tagged.csv > merged_counts_subset_noIDs.csv

#Merge both the subset and legacy gene count tables for further comparison
paste -d , merged_counts_legacy_tagged.csv merged_counts_subset_noIDs.csv > final_merged_counts.csv