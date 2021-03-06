#!/bin/bash
#!/bin/bash
#Script to re-format merged count tables to Gene Cluster Text file format (*.gct)
#Usage: bash reformatCounts_normalized_GCT.sh countsFolder sampleSet
#Usage Ex: bash reformatCounts_normalized_GCT.sh GeneCountsAnalyzed_countedCoordinate_htseqHisat2_run1_fullset R2

#Check for input argument of file name
if [ $# -eq 0 ]; then
   	echo "ERROR: No folder name(s) supplied... exiting"
   	exit 1
fi
#Retrieve statistics outputs absolute path
outputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Retrieve analysis inputs path
inputsPath=$(grep "statistics:" ../InputData/outputPaths.txt | tr -d " " | sed "s/statistics://g")
#Initialize values
numRows=0
numCols=0
#Retrieve input filename
inFile="$inputsPath"/GeneCounts_Stats/"$1"/geneCountStats_"$2"/stats_normalizedCounts.csv
#Change delimiter for input files
sed 's/,/\t/g' "$inFile" > tmpInFile.txt
sed -i 's/ /\t/g' tmpInFile.txt
#Retrieve number of rows
wc -l tmpInFile.txt > tmpNumRows.txt
numRows=$(cut -d ' ' -f1 tmpNumRows.txt)
#Retrieve number of samples
numCols=$(($(head -n1 tmpInFile.txt | awk '{print NF}')-1))
#Set output file name
outFile="$outputsPath"/GeneCounts_Stats/"$1"/geneCountStats_"$2"/stats_normalizedCounts_reformatted.gct
#Output headers for GCT formatting
echo "#1.2" > tmpHeader.gct
echo -e "$numRows \t $numCols" >> tmpHeader.gct
#Create temporary file with added empty second column for the 'description' field
cut -f1 tmpInFile.txt > tmpData1.gct
sed -e "s/$/\tNA/" -i tmpData1.gct
cut -f2- tmpInFile.txt > tmpData2.gct
paste tmpData1.gct tmpData2.gct > tmpData3.gct
sed -i "s/gene\tNA/Name\tDescription/g" tmpData3.gct
#Append header to reformatted counts table
cat tmpHeader.gct tmpData3.gct > "$outFile"
#Print a script completion confirmation message
echo "Gene counts file has been reformatted!"
#Clean up
rm tmp*.txt
rm tmp*.gct