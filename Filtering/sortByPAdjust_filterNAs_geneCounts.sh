#Script to filter and sort gene count files by adjusted p-value
#Usage: bash sortByPAdjust_filterNAs_geneCounts.sh geneCountsFileWithPAdjust
#Usage Ex: bash sortByPAdjust_filterNAs_geneCounts.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot_pAdjust.csv

#Retrieve inputs
inputCounts="$1"

#Set outputs file path
outputCounts=$(echo "$inputCounts" | sed "s/\.csv/\_cleaned\.csv/g")

#Create initial output file
cat $inputCounts > tmpData1.csv

#Remove all lines with NAs
sed -i ".bak" '/^NA,/d' tmpData1.csv
sed -i ".bak" '/,NA,/d' tmpData1.csv

#Add header to output file
head -1 tmpData1.csv > $outputCounts

#Sort the file by the adjusted p-value column
tail -n+2 tmpData1.csv | sort -k2 -t "," >> $outputCounts

#Clean up
rm tmp*