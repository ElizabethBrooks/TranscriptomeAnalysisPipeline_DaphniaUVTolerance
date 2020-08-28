#!/bin/bash
#Usage: bash prepareMatrices.sh countedGenesFile
#Usage Ex: bash prepareMatrices.sh geneCounts_merged_genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1.txt

#Prepare for analysis
dirFlag=0
runNum=1
#Trim file extension from input
inputTable=$(echo "$1" | sed "s/\.txt//g")
#Directory for outputs
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
outputCounts="$outputsPath"/GeneCountsAnalyzed
#Directory for gene count tables
inputCounts="$outputCounts"/"$1"

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$inputCounts" | tr -s '[:blank:]' ',' > "$outputCounts"/"$inputTable"_blankCleaned.csv
#Remove extra lines from the new gene count tables
egrep -v "unique|ambiguous|feature|aligned|aQual" "$outputCounts"/"$inputTable"_blankCleaned.csv > "$outputCounts"/"$inputTable"_cleaned.csv
#Clean up temporary files
rm "$outputCounts"/"$inputTable"_blankCleaned.csv

#Add postfix tags to each sample name in each file indicating
# the alignment method used (T for tophat and H for hisat2)
#Determine which analysis method was used
if [[ "$inputTable" == *"Hisat2"* ]]; then
	sed 's/Pool1/Pool1_H/g' "$outputCounts"/"$inputTable"_cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$outputCounts"/"$inputTable"_tagged.csv
else
	sed 's/Pool1/Pool1_T/g' "$outputCounts"/"$inputTable"_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$outputCounts"/"$inputTable"_tagged.csv
fi

#Transpose gene count tables for PCA and fix headers
csvtool transpose "$outputCounts"/"$inputTable"_cleaned.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/"$inputTable"_transposed.csv
csvtool transpose "$outputCounts"/"$inputTable"_tagged.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/"$inputTable"_tagged_transposed.csv

#Add column to transposed tables with treatment and alignment method before merging
# or add column to transposed tables with treatment
sed '1 s/$/,treatment/' "$outputCounts"/"$inputTable"_transposed.csv > "$outputCounts"/"$inputTable"_annotatedTreatment_transposed.csv
sed -i '/UV/ s/$/,UV/' "$outputCounts"/"$inputTable"_annotatedTreatment_transposed.csv
sed -i '/VIS/ s/$/,VIS/' "$outputCounts"/"$inputTable"_annotatedTreatment_transposed.csv

#Add column to transposed tables with geneotype
sed '1 s/$/,geneotype/' "$outputCounts"/"$inputTable"_transposed.csv > "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/Y05/ s/$/,Y05/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/Y023/ s/$/,Y023/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/E05/ s/$/,E05/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/R2/ s/$/,R2/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/PA/ s/$/,PA/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
sed -i '/Sierra/ s/$/,Sierra/' "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv
#Clean up
rm "$outputCounts"/"$inputTable"_tagged_transposed.csv

#Create merged count tables
#Transpose the annotated tables
csvtool transpose "$outputCounts"/"$inputTable"_annotatedTreatment_transposed.csv > "$outputCounts"/"$inputTable"_annotatedTreatment.csv
csvtool transpose "$outputCounts"/"$inputTable"_annotatedGeneotype_transposed.csv > "$outputCounts"/"$inputTable"_annotatedGeneotype.csv

#Print a script completion confirmation message
echo "Gene count matricies have been prepared!"