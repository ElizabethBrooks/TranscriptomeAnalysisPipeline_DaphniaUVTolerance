#Script to swap gene IDs with Uniprot IDs
#Usage: bash swapGeneToUniprot.sh countsFile
#Usage Ex: bash swapGeneToUniprot.sh ~/PfrenderLab/PA42_v4.1/geneCounts_cleaned_PA42_v4.1.csv
#Usage Ex: bash swapGeneToUniprot.sh ~/PfrenderLab/dMelUV/WGCNA_PA42_v4.1/normalizedCountsInter_PA42_v4.1.csv

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
cp $countPath $countFile

#Make tmp counts file without header
tail -n+2 $countFile > tmp.csv

#Loop over each gene ID in the counts file
while read -r line; do
	gTag=$(echo "$line" | cut -d"," -f1 | sed "s/\"//g")
	gTag=$gTag","
	uTag=$(grep "$gTag" $tagFile | cut -d"," -f2)
	#Swap gene with Uniprot ID, if there is a Uniprot ID
	if [[ "$uTag" != "." ]]; then
		sed -i ".bak" "s/$gTag/$uTag,/g" $countFile
	fi
done < tmp.csv

#Loop over each gene ID
#while read -r line; do
#	uTag=$(echo "$line" | cut -d"," -f2)
	#Swap gene with Uniprot ID, if there is a Uniprot ID
#	if [[ "$uTag" != "." ]]; then
#		gTag=$(echo "$line" | cut -d"," -f1)
#		sed -i ".bak" "s/$gTag,/$uTag,/g" $countFile
#	fi
#done < $tagFile

#Clean up
rm $countFile".bak"
rm tmp.csv