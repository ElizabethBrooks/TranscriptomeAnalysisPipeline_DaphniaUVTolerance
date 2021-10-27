#Script to add adjusted p-values to gene count files
# that have uniprot IDs added already
#Usage: bash addPAdjustToCounts.sh countsFile pValuesFile
#Usage Ex: bash addPAdjustToCounts.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot.csv /Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags_filtered.csv
#Usage Ex: bash addPAdjustToCounts.sh /Users/bamflappy/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1_uniprot.csv /Users/bamflappy/PfrenderLab/dMelUV/DEA_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags.csv

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve input gene count file path
countPath="$1"

#Set output gene count file path
countFile=$(echo "$countPath" | sed "s/\.csv/\_pAdjust\.csv/g")

#Set input tag file path
tagFile="$2"

#Separate first column of uniprot IDs
cut -f1 -d"," $countPath > tmpData1.csv

#Create temporary file with added empty column
cut -f2- -d"," $countPath > tmpData2.csv
sed -i '.bak' -e "s/^/NA,/" tmpData2.csv 

#Loop over each gene ID and add matching adjusted p-values
while read -r line; do
	gTag=$(echo "$line" | cut -d"," -f1)
	pTag=$(echo "$line" | cut -d"," -f6)
	sed -i ".bak" "s/NA,$gTag,/$pTag,$gTag,/g" tmpData2.csv
done < $tagFile

#Add column of uniprot IDs back to counts file
paste -d "," tmpData1.csv tmpData2.csv > $countFile

#Clean up
rm *.bak
rm tmp*.csv