#!/bin/bash
#Fullset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Subset gene counts were created using trimmomatic, hisat2, and cuffdiff
#Legacy gene counts were created using trimmomatic, tophat, and cuffdiff
#Retrieve gene count analysis outputs absolute path
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
#Move to outputs directory
cd "$outputsPath"
#Directory for gene count tables
prefixInputs="$outputsPath"/"GeneCounts_Merged/GeneCounts_Tables"
#Directory for outputs
prefixOutputs="$outputsPath"/"GeneCounts_Merged"

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$prefixInputs"/"$1" | tr -s '[:blank:]' ',' > "$prefixOutputs"/"$1"_blankCleaned.csv
cat "$prefixInputs"/merged_counts_legacy.txt | tr -s '[:blank:]' ',' > "$prefixOutputs"/merged_counts_legacy_cleaned.csv
#Remove extra lines from the new gene count tables
egrep -v "unique|ambiguous|feature|aligned|aQual" "$prefixOutputs"/"$1"_blankCleaned.csv > "$prefixOutputs"/"$1"_cleaned.csv
#Clean up temporary files
rm "$prefixOutputs"/"$1"_blankCleaned.csv

#Compare row tags to identify lines unique to either file,
# in order to ensure matching gene IDs before merging
#Retrieve the sorted row tags of gene IDs for each file
#cat "$prefixOutputs"/merged_counts_subset_cleaned.csv | cut -d, -f1 | sort -d > "$prefixOutputs"/merged_counts_subset_rowTags.csv
#cat "$prefixOutputs"/merged_counts_legacy_cleaned.csv | cut -d, -f1 | sort -d > "$prefixOutputs"/merged_counts_legacy_rowTags.csv
#Compare and supress output of row tags found in both files using the -3 flag
#comm -3 "$prefixOutputs"/merged_counts_subset_rowTags.csv "$prefixOutputs"/merged_counts_legacy_rowTags.csv > "$prefixOutputs"/merged_counts_uniqueRowTags.csv
#Clean up temporary files
#rm "$prefixOutputs"/merged_counts_subset_rowTags.csv
#rm "$prefixOutputs"/merged_counts_legacy_rowTags.csv

#Add postfix tags to each sample name in each file indicating
# the alignment method used (T for tophat and H for hisat2)
sed 's/Pool1/Pool1_H/g' "$prefixOutputs"/"$1"_cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$prefixOutputs"/"$1"_tagged.csv
sed 's/Pool1/Pool1_T/g' "$prefixOutputs"/merged_counts_legacy_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$prefixOutputs"/merged_counts_legacy_tagged.csv

#Transpose gene count tables for PCA and fix headers
#Fullset and Subset
csvtool transpose "$prefixOutputs"/"$1"_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$1"_transposed.csv
csvtool transpose "$prefixOutputs"/"$1"_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$1"_tagged_transposed.csv
#Legacy
csvtool transpose "$prefixOutputs"/merged_counts_legacy_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/merged_counts_legacy_transposed.csv
csvtool transpose "$prefixOutputs"/merged_counts_legacy_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv
#Clean up temporary files
rm "$prefixOutputs"/"$1"_cleaned.csv
rm "$prefixOutputs"/merged_counts_legacy_cleaned.csv

#Add column to transposed tables with treatment and alignment method before merging
# or add column to transposed tables with treatment
if [[ "$1" == *"subset"* ]]; then
	#Subset
	sed '1 s/$/,method/' "$prefixOutputs"/"$1"_tagged_transposed.csv > "$prefixOutputs"/merged_counts_subset_annotatedMethod_transposed.csv
	sed -i 's/$/,hisat2/' "$prefixOutputs"/"$1"_annotatedMethod_transposed.csv
	sed -i 's/\<method,hisat2\>/method/g' "$prefixOutputs"/"$1"_annotatedMethod_transposed.csv
	#Legacy
	sed '1 s/$/,method/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
	sed -i 's/$/,tophat/' "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
	sed -i 's/\<method,tophat\>/method/g' "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
	#Add column to transposed tables with treatment
	#Subset
	sed '1 s/$/,treatment/' "$prefixOutputs"/"$1"_tagged_transposed.csv > "$prefixOutputs"/merged_counts_subset_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv
	#Legacy
	sed '1 s/$/,treatment/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
else
	#Add column to transposed tables with treatment
	#Fullset
	sed '1 s/$/,treatment/' "$prefixOutputs"/"$1"_transposed.csv > "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv
fi

#Add column to transposed tables with geneotype
if [[ "$1" == *"subset"* ]]; then
	#Subset
	sed '1 s/$/,geneotype/' "$prefixOutputs"/"$1"_tagged_transposed.csv > "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	#Legacy
	sed '1 s/$/,geneotype/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
else
	#Fullset
	sed '1 s/$/,geneotype/' "$prefixOutputs"/"$1"_transposed.csv > "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/PA/ s/$/,PA/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
	sed -i '/Sierra/ s/$/,Sierra/' "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv
fi
#Clean up
rm "$prefixOutputs"/"$1"_tagged_transposed.csv
rm "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv

#Create merged count tables
#Transpose the annotated tables
if [[ "$1" == *"subset"* ]]; then
	#Subset
	csvtool transpose "$prefixOutputs"/"$1"_annotatedMethod_transposed.csv > "$prefixOutputs"/"$1"_annotatedMethod.csv
	csvtool transpose "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv > "$prefixOutputs"/"$1"_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv > "$prefixOutputs"/"$1"_annotatedGeneotype.csv
	#Legacy
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv
	#Clean up
	rm "$prefixOutputs"/"$1"_annotatedMethod_transposed.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
else
	#Fullset
	csvtool transpose "$prefixOutputs"/"$1"_annotatedTreatment_transposed.csv > "$prefixOutputs"/"$1"_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/"$1"_annotatedGeneotype_transposed.csv > "$prefixOutputs"/"$1"_annotatedGeneotype.csv
fi

if [[ "$1" == *"subset"* ]]; then
	#Remove the row tags with gene IDs from subset tables before merging
	#sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/merged_counts_subset_tagged.csv > "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$1"_annotatedMethod.csv > "$prefixOutputs"/"$1"_annotatedMethod_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$1"_annotatedTreatment.csv > "$prefixOutputs"/"$1"_annotatedTreatment_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$1"_annotatedGeneotype.csv > "$prefixOutputs"/"$1"_annotatedGeneotype_trimmed.csv
	#Clean up temporary files
	rm "$prefixOutputs"/"$1"_tagged.csv
	rm "$prefixOutputs"/"$1"_annotatedMethod.csv
	rm "$prefixOutputs"/"$1"_annotatedTreatment.csv
	rm "$prefixOutputs"/"$1"_annotatedGeneotype.csv

	#Merge both the subset and legacy gene count tables for further comparison
	#paste -d , "$prefixOutputs"/merged_counts_legacy_tagged.csv "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv > "$prefixOutputs"/"$1"_final_merged_counts_tagged.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv "$prefixOutputs"/"$1"_annotatedMethod_trimmed.csv > "$prefixOutputs"/"$1"_final_merged_counts_annotatedMethod.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv "$prefixOutputs"/"$1"_annotatedTreatment_trimmed.csv > "$prefixOutputs"/"$1"_final_merged_counts_annotatedTreatment.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv "$prefixOutputs"/"$1"_annotatedGeneotype_trimmed.csv > "$prefixOutputs"/"$1"_final_merged_counts_annotatedGeneotype.csv
	#Clean up temporary files
	rm "$prefixOutputs"/merged_counts_legacy_tagged.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv
	#rm "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv
	rm "$prefixOutputs"/"$1"_annotatedMethod_trimmed.csv
	rm "$prefixOutputs"/"$1"_annotatedTreatment_trimmed.csv
	rm "$prefixOutputs"/"$1"_annotatedGeneotype_trimmed.csv

	#Transpose merged count tables and fix headers
	#csvtool transpose "$prefixOutputs"/"$1"_final_merged_counts_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$1"_final_merged_counts_tagged_transposed.csv
	csvtool transpose "$prefixOutputs"/"$1"_final_merged_counts_annotatedMethod.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$1"_final_merged_counts_annotatedMethod_transposed.csv
	csvtool transpose "$prefixOutputs"/"$1"_final_merged_counts_annotatedTreatment.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$1"_final_merged_counts_annotatedTreatment_transposed.csv
	csvtool transpose "$prefixOutputs"/"$1"_final_merged_counts_annotatedGeneotype.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$1"_final_merged_counts_annotatedGeneotype_transposed.csv
	#Clean up temporary files
	#rm "$prefixOutputs"/"$1"_final_merged_counts_tagged.csv
	rm "$prefixOutputs"/"$1"_final_merged_counts_annotatedMethod.csv
	rm "$prefixOutputs"/"$1"_final_merged_counts_annotatedTreatment.csv
	rm "$prefixOutputs"/"$1"_final_merged_counts_annotatedGeneotype.csv
fi