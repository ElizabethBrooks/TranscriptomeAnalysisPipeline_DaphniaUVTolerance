#!/bin/bash
# usage: bash formatCounts.sh
# usage Ex: bash formatCounts.sh 

#Trim file extension from input
inputCounts=$(grep "geneCounts:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneCounts://g")

#Directory for outputs
#outputsPath=$(dirname "$inputCounts")
outputsPath=$(grep "cleanedGeneCounts:" ../InputData/outputPaths.txt | tr -d " " | sed "s/cleanedGeneCounts://g")

# set outputs path
#outputCounts=$outputsPath"/Formatted"
outputCounts=$outputsPath

#Make outputs directory
mkdir $outputCounts

#Prepare new gene count tables for comparison
#Remove excess white space in gene count tables,
# particularly a sample name in the header
cat "$inputCounts" | tr -s '[:blank:]' ',' > "$outputCounts"/pre_cleaned.csv
#Remove extra lines from the new gene count tables
#egrep -v "unique|ambiguous|feature|aligned|aQual" "$outputCounts"/blankCleaned.csv > "$outputCounts"/pre_cleaned.csv

#Add postfix tags to each sample name in each file indicating
# the alignment method used (T for tophat and H for hisat2)
#Determine which analysis method was used
#if [[ "$inputCounts" == *"Hisat2"* ]]; then
sed 's/Pool1/Pool1_H/g' "$outputCounts"/pre_cleaned.csv |sed 's/Pool2/Pool2_H/g' | sed 's/Pool3/Pool3_H/g' > "$outputCounts"/tagged.csv
#elif [[ "$inputCounts" == *"Tophat2"*  ]]; then
	#statements
#	sed 's/Pool1/Pool1_T/g' "$outputCounts"/pre_cleaned.csv |sed 's/Pool2/Pool2_T/g' | sed 's/Pool3/Pool3_T/g' > "$outputCounts"/tagged.csv
#else
#	echo "A valid alignment method was not detected... exiting!"
#	exit 1
#fi

#Transpose gene count tables and fix headers
csvtk transpose "$outputCounts"/pre_cleaned.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/transposed.csv
csvtk transpose "$outputCounts"/tagged.csv | sed 's/\<gene\>/sample/g' > "$outputCounts"/tagged_transposed.csv

#Add column to transposed tables with treatment and alignment method before merging
# or add column to transposed tables with treatment
sed '1 s/$/,treatment/' "$outputCounts"/transposed.csv > "$outputCounts"/annotatedTreatment_transposed.csv
sed -i.bu '/UV/ s/$/,UV/' "$outputCounts"/annotatedTreatment_transposed.csv
sed -i.bu '/VIS/ s/$/,VIS/' "$outputCounts"/annotatedTreatment_transposed.csv

#Add column to transposed tables with geneotype
sed '1 s/$/,geneotype/' "$outputCounts"/transposed.csv > "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/Y05/ s/$/,Y05/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/Y023/ s/$/,Y023/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/E05/ s/$/,E05/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/R2/ s/$/,R2/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/PA/ s/$/,PA/' "$outputCounts"/annotatedGeneotype_transposed.csv
sed -i.bu '/Sierra/ s/$/,Sierra/' "$outputCounts"/annotatedGeneotype_transposed.csv
#Clean up
rm "$outputCounts"/tagged_transposed.csv
rm "$outputCounts"/*.csv.bu

#Create merged count tables
#Transpose the annotated tables
csvtk transpose "$outputCounts"/annotatedTreatment_transposed.csv > "$outputCounts"/annotatedTreatment.csv
csvtk transpose "$outputCounts"/annotatedGeneotype_transposed.csv > "$outputCounts"/annotatedGeneotype.csv

# remove gene: tags
cat "$outputCounts"/pre_cleaned.csv | sed "s/gene\://g" > "$outputCounts"/cleaned.csv

#Print a script completion confirmation message
echo "Gene count matricies have been prepared!"