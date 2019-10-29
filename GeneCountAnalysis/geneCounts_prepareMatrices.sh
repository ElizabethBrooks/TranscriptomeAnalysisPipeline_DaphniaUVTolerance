#!/bin/bash
#Fullset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Subset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Legacy gene counts were created using trimmomatic, tophat, and cuffdiff

#Directory for gene count tables
prefixInputs="../../GeneCounts_Merged/GeneCounts_Tables/"
#Directory for outputs
prefixOutputs="../../GeneCounts_Merged/"

#Merge counts from fullset of new gene count table
#Add directory paths to each file name in the guide file
#sed -i -e 's|^|stats_hisat2EdgeR_test1/stats_hisat2EdgeR_counts_test1/|' "$prefixOutputs"mergeGuideFile_fullset.txt
#python "$prefixOutputs"merge_tables.py "$prefixOutputs"mergeGuideFile_fullset.txt
#Rename the output merged counts file
#mv "$prefixOutputs"merged_counts.txt "$prefixOutputs"merged_counts_fullset.txt

#Merge counts from subset of new gene count table matching samples
# in the previously created "legacy" gene count table
#Add directory paths to each file name in the guide file
#sed -i -e 's|^|stats_hisat2EdgeR_test1/stats_hisat2EdgeR_counts_test1/|' "$prefixOutputs"mergeGuideFile_subset.txt
#python "$prefixOutputs"merge_tables.py "$prefixOutputs"mergeGuideFile_subset.txt
#Rename the output merged counts file
#mv "$prefixOutputs"merged_counts.txt "$prefixOutputs"merged_counts_subset.txt

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$prefixInputs"merged_counts_fullset.txt | tr -s '[:blank:]' ',' > "$prefixOutputs"merged_counts_fullset_blankCleaned.csv
cat "$prefixInputs"merged_counts_subset.txt | tr -s '[:blank:]' ',' > "$prefixOutputs"merged_counts_subset_blankCleaned.csv
cat "$prefixInputs"merged_counts_legacy.txt | tr -s '[:blank:]' ',' > "$prefixOutputs"merged_counts_legacy_cleaned.csv
#Remove extra lines from the new gene count tables
egrep -v "unique|ambiguous|feature|aligned|aQual" "$prefixOutputs"merged_counts_fullset_blankCleaned.csv > "$prefixOutputs"merged_counts_fullset_cleaned.csv
egrep -v "unique|ambiguous|feature|aligned|aQual" "$prefixOutputs"merged_counts_subset_blankCleaned.csv > "$prefixOutputs"merged_counts_subset_cleaned.csv
#Clean up temporary files
rm "$prefixOutputs"merged_counts_fullset_blankCleaned.csv
rm "$prefixOutputs"merged_counts_subset_blankCleaned.csv

#Compare row tags to identify lines unique to either file,
# in order to ensure matching gene IDs before merging
#Retrieve the sorted row tags of gene IDs for each file
#cat "$prefixOutputs"merged_counts_subset_cleaned.csv | cut -d, -f1 | sort -d > "$prefixOutputs"merged_counts_subset_rowTags.csv
#cat "$prefixOutputs"merged_counts_legacy_cleaned.csv | cut -d, -f1 | sort -d > "$prefixOutputs"merged_counts_legacy_rowTags.csv
#Compare and supress output of row tags found in both files using the -3 flag
#comm -3 "$prefixOutputs"merged_counts_subset_rowTags.csv "$prefixOutputs"merged_counts_legacy_rowTags.csv > "$prefixOutputs"merged_counts_uniqueRowTags.csv
#Clean up temporary files
#rm "$prefixOutputs"merged_counts_subset_rowTags.csv
#rm "$prefixOutputs"merged_counts_legacy_rowTags.csv

#Add postfix tags to each sample name in each file indicating
# the alignment method used (T for tophat and H for hisat2)
sed 's/Pool1/Pool1_H/g' "$prefixOutputs"merged_counts_subset_cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$prefixOutputs"merged_counts_subset_tagged.csv
sed 's/Pool1/Pool1_T/g' "$prefixOutputs"merged_counts_legacy_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$prefixOutputs"merged_counts_legacy_tagged.csv

#Transpose gene count tables for PCA and fix headers
#Fullset
csvtool transpose "$prefixOutputs"merged_counts_fullset_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"merged_counts_fullset_transposed.csv
#Subset
csvtool transpose "$prefixOutputs"merged_counts_subset_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"merged_counts_subset_transposed.csv
csvtool transpose "$prefixOutputs"merged_counts_subset_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"merged_counts_subset_tagged_transposed.csv
#Legacy
csvtool transpose "$prefixOutputs"merged_counts_legacy_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"merged_counts_legacy_transposed.csv
csvtool transpose "$prefixOutputs"merged_counts_legacy_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"merged_counts_legacy_tagged_transposed.csv
#Clean up temporary files
#rm "$prefixOutputs"merged_counts_fullset_cleaned.csv
#rm "$prefixOutputs"merged_counts_subset_cleaned.csv
#rm "$prefixOutputs"merged_counts_legacy_cleaned.csv

#Add column to transposed tables with treatment and alignment method
#Fullset
sed '1 s/$/,treatment/' "$prefixOutputs"merged_counts_fullset_transposed.csv > "$prefixOutputs"merged_counts_fullset_annotatedTreatment_transposed.csv
sed -i '/UV/ s/$/,UV/' "$prefixOutputs"merged_counts_fullset_annotatedTreatment_transposed.csv
sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"merged_counts_fullset_annotatedTreatment_transposed.csv
#Subset
sed '1 s/$/,treatment/' "$prefixOutputs"merged_counts_subset_tagged_transposed.csv > "$prefixOutputs"merged_counts_subset_annotatedTreatment_transposed.csv
sed -i '/UV/ s/$/,UV/' "$prefixOutputs"merged_counts_subset_annotatedTreatment_transposed.csv
sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"merged_counts_subset_annotatedTreatment_transposed.csv
sed '1 s/$/,method/' "$prefixOutputs"merged_counts_subset_tagged_transposed.csv > "$prefixOutputs"merged_counts_subset_annotatedMethod_transposed.csv
sed -i 's/$/,hisat2/' "$prefixOutputs"merged_counts_subset_annotatedMethod_transposed.csv
sed -i 's/\<method,hisat2\>/method/g' "$prefixOutputs"merged_counts_subset_annotatedMethod_transposed.csv
#Legacy
sed '1 s/$/,treatment/' "$prefixOutputs"merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"merged_counts_legacy_annotatedTreatment_transposed.csv
sed -i '/UV/ s/$/,UV/' "$prefixOutputs"merged_counts_legacy_annotatedTreatment_transposed.csv
sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"merged_counts_legacy_annotatedTreatment_transposed.csv
sed '1 s/$/,method/' "$prefixOutputs"merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"merged_counts_legacy_annotatedMethod_transposed.csv
sed -i 's/$/,tophat/' "$prefixOutputs"merged_counts_legacy_annotatedMethod_transposed.csv
sed -i 's/\<method,tophat\>/method/g' "$prefixOutputs"merged_counts_legacy_annotatedMethod_transposed.csv
#Clean up
rm "$prefixOutputs"merged_counts_subset_tagged_transposed.csv
rm "$prefixOutputs"merged_counts_legacy_tagged_transposed.csv

#Create merged count tables
#Transpose the annotated tables
#Fullset
csvtool transpose "$prefixOutputs"merged_counts_fullset_annotatedTreatment_transposed.csv > "$prefixOutputs"merged_counts_fullset_annotatedTreatment.csv
#Subset
csvtool transpose "$prefixOutputs"merged_counts_subset_annotatedTreatment_transposed.csv > "$prefixOutputs"merged_counts_subset_annotatedTreatment.csv
csvtool transpose "$prefixOutputs"merged_counts_subset_annotatedMethod_transposed.csv > "$prefixOutputs"merged_counts_subset_annotatedMethod.csv
#Legacy
csvtool transpose "$prefixOutputs"merged_counts_legacy_annotatedTreatment_transposed.csv > "$prefixOutputs"merged_counts_legacy_annotatedTreatment.csv
csvtool transpose "$prefixOutputs"merged_counts_legacy_annotatedMethod_transposed.csv > "$prefixOutputs"merged_counts_legacy_annotatedMethod.csv
#Clean up
rm "$prefixOutputs"merged_counts_subset_annotatedMethod_transposed.csv
rm "$prefixOutputs"merged_counts_legacy_annotatedMethod_transposed.csv

#Remove the row tags with gene IDs from subset tables before merging
#sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"merged_counts_subset_tagged.csv > "$prefixOutputs"merged_counts_subset_tagged_trimmed.csv
sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"merged_counts_subset_annotatedTreatment.csv > "$prefixOutputs"merged_counts_subset_annotatedTreatment_trimmed.csv
sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"merged_counts_subset_annotatedMethod.csv > "$prefixOutputs"merged_counts_subset_annotatedMethod_trimmed.csv
#Clean up temporary files
rm "$prefixOutputs"merged_counts_subset_tagged.csv
rm "$prefixOutputs"merged_counts_subset_annotatedTreatment.csv
rm "$prefixOutputs"merged_counts_subset_annotatedMethod.csv

#Merge both the subset and legacy gene count tables for further comparison
#paste -d , "$prefixOutputs"merged_counts_legacy_tagged.csv "$prefixOutputs"merged_counts_subset_tagged_trimmed.csv > "$prefixOutputs"final_merged_counts_tagged.csv
paste -d , "$prefixOutputs"merged_counts_legacy_annotatedTreatment.csv "$prefixOutputs"merged_counts_subset_annotatedTreatment_trimmed.csv > "$prefixOutputs"final_merged_counts_annotatedTreatment.csv
paste -d , "$prefixOutputs"merged_counts_legacy_annotatedMethod.csv "$prefixOutputs"merged_counts_subset_annotatedMethod_trimmed.csv > "$prefixOutputs"final_merged_counts_annotatedMethod.csv
#Clean up temporary files
rm "$prefixOutputs"merged_counts_legacy_tagged.csv
rm "$prefixOutputs"merged_counts_legacy_annotatedTreatment.csv
rm "$prefixOutputs"merged_counts_legacy_annotatedMethod.csv
#rm "$prefixOutputs"merged_counts_subset_tagged_trimmed.csv
rm "$prefixOutputs"merged_counts_subset_annotatedTreatment_trimmed.csv
rm "$prefixOutputs"merged_counts_subset_annotatedMethod_trimmed.csv

#Transpose merged count tables and fix headers
#csvtool transpose "$prefixOutputs"final_merged_counts_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"final_merged_counts_tagged_transposed.csv
csvtool transpose "$prefixOutputs"final_merged_counts_annotatedTreatment.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"final_merged_counts_annotatedTreatment_transposed.csv
csvtool transpose "$prefixOutputs"final_merged_counts_annotatedMethod.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"final_merged_counts_annotatedMethod_transposed.csv
#Clean up temporary files
#rm "$prefixOutputs"final_merged_counts_tagged.csv
rm "$prefixOutputs"final_merged_counts_annotatedTreatment.csv
rm "$prefixOutputs"final_merged_counts_annotatedMethod.csv