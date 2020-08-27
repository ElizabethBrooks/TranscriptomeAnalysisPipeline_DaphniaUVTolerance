#!/bin/bash
#Usage: bash prepareMatrices.sh countedGenesFile
#Usage Ex: bash prepareMatrices.sh geneCounts_merged_counted_htseqTophat2_run1_fullset.txt
#Fullset gene counts were created using trimmomatic, hisat2 or tophat2, and cuffdiff
#Subset gene counts were created using trimmomatic, hisat2 or tophat2, and cuffdiff
#Legacy gene counts were created using trimmomatic, tophat2, and cuffdiff
#Prepare for analysis
dirFlag=0
runNum=1
#Trim file extension from input
inputTable=$(echo "$1" | sed "s/\.txt//g")
#Directory for outputs
outputsPath=$(grep "geneTableAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneTableAnalysis://g")
prefixOutputs="$outputsPath"/GeneCounts_Formatted
#Directory for gene count tables
prefixInputs="$outputsPath"/GeneCounts_Formatted/"GeneCounts_Tables"
#Check which set of data is used
if [[ "$inputTable" == *"subset"* ]]; then
	analysisSet="subset"
else
	analysisSet="fullset"
fi
#Make a new directory for each analysis run
while [ $dirFlag -eq 0 ]; do
	prefixOutputs="$prefixOutputs"/GeneCountAnalysis_"$analysisSet"_run"$runNum"
	mkdir $prefixOutputs
	#Check if the folder already exists
	if [ $? -ne 0 ]; then
		#Increment the folder name
		let runNum+=1
	else
		#Indicate that the folder was successfully made
		dirFlag=1
		echo "Creating folder for run $runNum of gene count analysis of $1..."
	fi
done

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$prefixInputs"/"$inputTable".txt | tr -s '[:blank:]' ',' > "$prefixOutputs"/"$inputTable"_blankCleaned.csv
cat "$prefixInputs"/merged_counts_legacy.txt | tr -s '[:blank:]' ',' > "$prefixOutputs"/merged_counts_legacy_cleaned.csv
#Remove extra lines from the new gene count tables
egrep -v "unique|ambiguous|feature|aligned|aQual" "$prefixOutputs"/"$inputTable"_blankCleaned.csv > "$prefixOutputs"/"$inputTable"_cleaned.csv
#Clean up temporary files
rm "$prefixOutputs"/"$inputTable"_blankCleaned.csv

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
#Determine which analysis method was used
if [[ "$inputTable" == *"Hisat2"* ]]; then
	sed 's/Pool1/Pool1_H/g' "$prefixOutputs"/"$inputTable"_cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$prefixOutputs"/"$inputTable"_tagged.csv
else
	#Subset
	sed 's/Pool1/Pool1_T/g' "$prefixOutputs"/"$inputTable"_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$prefixOutputs"/"$inputTable"_tagged.csv
	#Legacy
	sed 's/Pool1/Pool1_T/g' "$prefixOutputs"/merged_counts_legacy_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$prefixOutputs"/merged_counts_legacy_tagged.csv
fi

#Transpose gene count tables for PCA and fix headers
#Fullset and Subset
csvtool transpose "$prefixOutputs"/"$inputTable"_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$inputTable"_transposed.csv
csvtool transpose "$prefixOutputs"/"$inputTable"_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$inputTable"_tagged_transposed.csv
if [[ "$inputTable" == *"subset"* ]]; then
	#Legacy
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_cleaned.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/merged_counts_legacy_transposed.csv
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv
	#Clean up temporary files
	#rm "$prefixOutputs"/merged_counts_legacy_cleaned.csv
fi
#Clean up temporary files
#rm "$prefixOutputs"/"$inputTable"_cleaned.csv

#Add column to transposed tables with treatment and alignment method before merging
# or add column to transposed tables with treatment
#Determine which set of data is being analyzed
if [[ "$inputTable" == *"subset"* ]]; then
	#Determine which analysis method was used
	if [[ "$inputTable" == *"Hisat2"* ]]; then
		sed '1 s/$/,method/' "$prefixOutputs"/"$inputTable"_tagged_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
		sed -i 's/$/,hisat2/' "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
		sed -i 's/\<method,hisat2\>/method/g' "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
	else
		#Subset
		sed '1 s/$/,method/' "$prefixOutputs"/"$inputTable"_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
		sed -i 's/$/,tophat/' "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
		sed -i 's/\<method,tophat\>/method/g' "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
		#Legacy
		sed '1 s/$/,method/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
		sed -i 's/$/,tophat/' "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
		sed -i 's/\<method,tophat\>/method/g' "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
	fi
	#Add column to transposed tables with treatment
	#Subset
	sed '1 s/$/,treatment/' "$prefixOutputs"/"$inputTable"_tagged_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
	#Legacy
	sed '1 s/$/,treatment/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv
else
	#Add column to transposed tables with treatment
	#Fullset
	sed '1 s/$/,treatment/' "$prefixOutputs"/"$inputTable"_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
	sed -i '/UV/ s/$/,UV/' "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
	sed -i '/VIS/ s/$/,VIS/' "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv
fi

#Add column to transposed tables with geneotype
if [[ "$inputTable" == *"subset"* ]]; then
	#Subset
	sed '1 s/$/,geneotype/' "$prefixOutputs"/"$inputTable"_tagged_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	#Legacy
	sed '1 s/$/,geneotype/' "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv
else
	#Fullset
	sed '1 s/$/,geneotype/' "$prefixOutputs"/"$inputTable"_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/Y05/ s/$/,Y05/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/Y023/ s/$/,Y023/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/E05/ s/$/,E05/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/R2/ s/$/,R2/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/PA/ s/$/,PA/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	sed -i '/Sierra/ s/$/,Sierra/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv
	#Clean up
	rm "$prefixOutputs"/merged_counts_legacy_tagged_transposed.csv
fi
#Clean up
rm "$prefixOutputs"/"$inputTable"_tagged_transposed.csv

#Create merged count tables
#Transpose the annotated tables
if [[ "$inputTable" == *"subset"* ]]; then
	#Subset
	csvtool transpose "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedMethod.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedGeneotype.csv
	#Legacy
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedTreatment_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype_transposed.csv > "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv
	#Clean up
	rm "$prefixOutputs"/"$inputTable"_annotatedMethod_transposed.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedMethod_transposed.csv
else
	#Fullset
	csvtool transpose "$prefixOutputs"/"$inputTable"_annotatedTreatment_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedTreatment.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_annotatedGeneotype_transposed.csv > "$prefixOutputs"/"$inputTable"_annotatedGeneotype.csv
fi

if [[ "$inputTable" == *"subset"* ]]; then
	#Remove the row tags with gene IDs from subset tables before merging
	#sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/merged_counts_subset_tagged.csv > "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$inputTable"_annotatedMethod.csv > "$prefixOutputs"/"$inputTable"_annotatedMethod_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$inputTable"_annotatedTreatment.csv > "$prefixOutputs"/"$inputTable"_annotatedTreatment_trimmed.csv
	sed 's/\([^,]*\),\(.*\)/\2/' "$prefixOutputs"/"$inputTable"_annotatedGeneotype.csv > "$prefixOutputs"/"$inputTable"_annotatedGeneotype_trimmed.csv
	#Clean up temporary files
	rm "$prefixOutputs"/"$inputTable"_tagged.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedMethod.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedTreatment.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedGeneotype.csv

	#Merge both the subset and legacy gene count tables for further comparison
	#paste -d , "$prefixOutputs"/merged_counts_legacy_tagged.csv "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv > "$prefixOutputs"/"$inputTable"_final_merged_counts_tagged.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv "$prefixOutputs"/"$inputTable"_annotatedMethod_trimmed.csv > "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedMethod.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv "$prefixOutputs"/"$inputTable"_annotatedTreatment_trimmed.csv > "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedTreatment.csv
	paste -d , "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv "$prefixOutputs"/"$inputTable"_annotatedGeneotype_trimmed.csv > "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedGeneotype.csv
	#Clean up temporary files
	rm "$prefixOutputs"/merged_counts_legacy_tagged.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedMethod.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedTreatment.csv
	rm "$prefixOutputs"/merged_counts_legacy_annotatedGeneotype.csv
	#rm "$prefixOutputs"/merged_counts_subset_tagged_trimmed.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedMethod_trimmed.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedTreatment_trimmed.csv
	rm "$prefixOutputs"/"$inputTable"_annotatedGeneotype_trimmed.csv

	#Transpose merged count tables and fix headers
	#csvtool transpose "$prefixOutputs"/"$inputTable"_final_merged_counts_tagged.csv | sed 's/\<gene\>/sample/g' > "$prefixOutputs"/"$inputTable"_final_merged_counts_tagged_transposed.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedMethod.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedMethod_transposed.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedTreatment.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedTreatment_transposed.csv
	csvtool transpose "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedGeneotype.csv | sed 's/\<gene\>/sample/g'> "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedGeneotype_transposed.csv
	#Clean up temporary files
	#rm "$prefixOutputs"/"$inputTable"_final_merged_counts_tagged.csv
	rm "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedMethod.csv
	rm "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedTreatment.csv
	rm "$prefixOutputs"/"$inputTable"_final_merged_counts_annotatedGeneotype.csv
fi
#Print a script completion confirmation message
echo "Gene count matricies have been prepared!"