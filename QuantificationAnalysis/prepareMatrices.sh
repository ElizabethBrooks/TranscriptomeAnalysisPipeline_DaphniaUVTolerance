#!/bin/bash
#Usage: bash prepareMatrices.sh countedGenesFile
#Usage Ex: bash prepareMatrices.sh geneCounts_merged_genome_sortedName_samtoolsHisat2_run1_counted_htseq_run1.txt

#Prepare for analysis
dirFlag=0
runNum=1
#Trim file extension from input
inputTable=$(echo "$1" | sed "s/\.txt//g" | sed "s/geneCounts\_merged\_//g")
#Directory for outputs
outputsPath=$(grep "geneCountAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/geneCountAnalysis://g")
outputCounts="$outputsPath"/GeneCountsAnalyzed/"$inputTable"
#Directory for gene count tables
inputCounts="$outputsPath"/GeneCountsAnalyzed/"$1"

#Make output directory
mkdir "$outputCounts"

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$inputCounts" | tr -s '[:blank:]' ',' > "$outputCounts"/blankCleaned.csv
#Remove extra lines from the new gene count tables
egrep -v "unique|ambiguous|feature|aligned|aQual" "$outputCounts"/blankCleaned.csv > "$outputCounts"/cleaned.csv
#Clean up temporary files
rm "$outputCounts"/blankCleaned.csv

#Add postfix tags to each sample name in each file indicating
# the alignment method used (T for tophat and H for hisat2)
#Determine which analysis method was used
if [[ "$inputTable" == *"Hisat2"* ]]; then
	sed 's/Pool1/Pool1_H/g' "$outputCounts"/cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$outputCounts"/tagged.csv
else
	sed 's/Pool1/Pool1_T/g' "$outputCounts"/cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$outputCounts"/tagged.csv
fi

#Transpose gene count tables for PCA and fix headers
csvtool transpose "$outputCounts"/cleaned.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/transposed.csv
csvtool transpose "$outputCounts"/tagged.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/tagged_transposed.csv

#Add column to transposed tables with treatment and alignment method before merging
# or add column to transposed tables with treatment
sed '1 s/$/,treatment/' "$outputCounts"/transposed.csv > "$outputCounts"/annotatedTreatment_transposed.csv
sed -i '/UV/ s/$/,UV/' "$outputCounts"/annotatedTreatment_transposed.csv
sed -i '/VIS/ s/$/,VIS/' "$outputCounts"/annotatedTreatment_transposed.csv

#Add column to transposed tables with geneotype
sed '1 s/$/,geneotype/' "$outputCounts"/transposed.csv > "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/Y05/ s/$/,Y05/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/Y023/ s/$/,Y023/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/E05/ s/$/,E05/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/R2/ s/$/,R2/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/PA/ s/$/,PA/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i '/Sierra/ s/$/,Sierra/' "$outputCounts"/annotatedGeneotype_transposed.csv
#Clean up
rm "$outputCounts"/tagged_transposed.csv

#Create merged count tables
#Transpose the annotated tables
csvtool transpose "$outputCounts"/annotatedTreatment_transposed.csv > "$outputCounts"/annotatedTreatment.csv
csvtool transpose "$outputCounts"/annotatedGeneotype_transposed.csv > "$outputCounts"/annotatedGeneotype.csv

#Print a script completion confirmation message
echo "Gene count matricies have been prepared!"