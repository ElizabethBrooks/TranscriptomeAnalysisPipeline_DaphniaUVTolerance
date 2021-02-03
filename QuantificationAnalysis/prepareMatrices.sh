#!/bin/bash
#Usage: bash prepareMatrices.sh countedGenesFilePath
#Usage: bash prepareMatrices.sh /home/mae/Documents/RNASeq_Workshop_ND/geneCounts_merged_genome_sortedName_samtoolsHisat2_run2_counted_htseq_run1.txt

#Prepare for analysis
dirFlag=0
runNum=1
#Trim file extension from input
inputCounts="$1"
#Directory for outputs
outputsPath=$(dirname "$inputCounts")
outputCounts="$outputsPath"

#Make output directory
inputTag=$(basename "$inputCounts" | sed 's/\.txt//g')
outputCounts="$outputCounts"/"$inputTag"
mkdir "$outputCounts"
#Check if the folder already exists
if [ $? -ne 0 ]; then
	echo "The $outputCounts directory already exsists... please remove before proceeding."
	exit 1
fi

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
if [[ "$inputCounts" == *"Hisat2"* ]]; then
	sed 's/Pool1/Pool1_H/g' "$outputCounts"/cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$outputCounts"/tagged.csv
elif [[ "$inputCounts" == *"Tophat2"*  ]]; then
	#statements
	sed 's/Pool1/Pool1_T/g' "$outputCounts"/cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$outputCounts"/tagged.csv
else
	echo "A valid alignment method was not detected... exiting!"
	exit 1
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