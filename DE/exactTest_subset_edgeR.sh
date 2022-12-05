#!/bin/bash
#Script to run Rscripts that perform DE analysis of gene count tables
#Usage: bash exactTest_subset_edgeR.sh sample FDR
#Usage Ex: bash exactTest_subset_edgeR.sh R2 0.10

#Check for input arguments of folder names
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi

#Retrieve gene ontology data path
#ontologyPath=$(grep "geneOntology:" ../InputData/inputPaths.txt | tr -d " " | sed "s/geneOntology://g")

#Retrieve analysis outputs path
outputsPath=$(grep "DEAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/DEAnalysis://g")

#Retrieve analysis inputs path
#inFile=$(grep "cleanedGeneCounts:" ../InputData/outputPaths_"$2".txt | tr -d " " | sed "s/cleanedGeneCounts://g")
inFile=$(grep "cleanedGeneCounts:" ../InputData/outputPaths.txt | tr -d " " | sed "s/cleanedGeneCounts://g")
inFile=$inFile"/cleaned.csv"

#Set FDR cut off
fdrCut=$2

#Create directory for output files
outDir="$outputsPath"/exactTest"$1"Analysis_FDR"$fdrCut"
mkdir $outDir

#Convert TXT formatted counts to CSV
#sed -e 's/\s\+/,/g' "$inFile" > "$newFile"
#cat "$inFile" > "$newFile"

#Retrieve selected sample beginning column number
#head -1 "$inFile" | tr "\t" "\n" | grep -n "$1"
colNumStart=$(($(head -1 "$inFile" | tr "," "\n" | grep -n "$1" | head -1 | cut -d ':' -f1)-1))
colNumEnd=$(($colNumStart+5))

#Perform DE analysis using edgeR and output analysis results to a txt file
Rscript exactTest_edgeR.r "$inFile" $colNumStart $colNumEnd $fdrCut > "$outDir"/analysisResults.txt

#Fix formatting and headers for the normalized counts table and exact test stats
#sed -i 's/"//g' stats_normalizedCounts.csv
#sed -i 's/"//g' stats_exactTest.csv
#sed -i -e 's/^$2_VIS_Pool1/gene,$2_VIS_Pool1/' stats_normalizedCounts.csv
#Add gene tags, mRNA tags for assemblies
#if [[ "$1" == genome* ]]; then
#	sed -i -e 's/^logFC/gene,logFC/' stats_exactTest.csv
#else
#	sed -i -e 's/^logFC/mRNA,logFC/' stats_exactTest.csv
#fi

#Move produced tables
for f in *.csv; do
	file="$f"
	sed -i 's/"//g' "$file"
	#Fix header
	tail -n+2 "$file" > tmpTail.csv
	head -1 "$file" | sed -e 's/^/gene,/' > tmpHeader.csv
	#Update table
	cat tmpHeader.csv > "$file"
	cat tmpTail.csv >> "$file"
	rm tmp*.csv
	#Move updated table
	mv "$file" "$outDir"
done
#Move produced plots
for p in *.jpg; do mv $p "$outDir"; done

#Move to current outputs folder
#cd "$outDir"
#Make table of GO data for the top tags from exact tests
#head -11 stats_exactTest.csv > topGenesStats_exactTest.csv
#sed -i 's/^logFC/gene,logFC/g' topGenesStats_exactTest.csv
#cut -f1 -d ',' topGenesStats_exactTest.csv > tmp.csv
#head -1 "$ontologyPath" > topGenesGO_exactTest.csv
#while IFS= read -r line; do grep $line "$ontologyPath" >> topGenesGO_exactTest.csv; done < "tmp.csv"
#Clean up
#rm tmp.csv