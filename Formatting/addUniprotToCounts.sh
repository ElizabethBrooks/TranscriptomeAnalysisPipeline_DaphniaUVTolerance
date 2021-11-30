#Script to swap gene IDs with Uniprot IDs
#Usage: bash addUniprotToCounts.sh countsFile
#Usage Ex: bash addUniprotToCounts.sh ~/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv
#Usage Ex: bash addUniprotToCounts.sh ~/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/filteredCountFiles/normalizedCountsInter_PA42_v4.1.csv
#Usage Ex: bash addUniprotToCounts.sh ~/PfrenderLab/PA42_v4.1/PA42_v4.1_normalizedCountsOlympics.csv

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve input gene count file path
countPath="$1"

#Set output gene count file path
countFile=$(echo "$countPath" | sed "s/\.csv/\_uniprot\.csv/g")

#Set input tag file path
tagFile="/Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_geneUniprotMap_cleaned.csv"

#Copy the gene count file to the results file
cat $countPath | sed 's/"",/gene,/g' | sed "s/\"//g" | sed 's/^Y05_VIS_Pool1,/gene,Y05_VIS_Pool1,/g' > $countFile

#Add initial column of NAs for Uniprot IDs
sed -i '.bak' -e "s/^/NoKnownID,/" $countFile

#Loop over each gene ID
while read -r line; do
	gTag=$(echo "$line" | cut -d"," -f1)
	uTag=$(echo "$line" | cut -d"," -f2)
	#Swap gene with Uniprot ID, if there is a Uniprot ID
	#Add the gene IDs back in as the description column
	if [[ "$uTag" != "." ]]; then
		sed -i ".bak" "s/NoKnownID,$gTag,/$uTag,$gTag,/g" $countFile
	else
		sed -i ".bak" "s/NoKnownID,$gTag,/NoKnownID_$gTag,$gTag,/g" $countFile
	fi
done < $tagFile

#Clean up
rm $countFile".bak"