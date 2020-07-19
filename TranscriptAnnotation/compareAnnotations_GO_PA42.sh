#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations_GO_PA42.sh annotationFile1 annotationFile2
#Usage ex: bash compareAnnotations_GO_PA42.sh PA42_proteins/GO.out.txt PA42_cds/GO.out.txt

echo "Beginning GO comparisons..."

#Retrieve ID fields
tail -n +2 $1 | cut -f 1,3 > tmp1_fields.csv
tail -n +2 $2 | cut -f 1,3 > tmp2_fields.csv

#Change white spaces to commas
# and change mRNA tags to gene
sed -e 's/\s\+/,/g' tmp1_fields.csv | sed 's/mRNA/gene/g' | sed 's/\.p1//g' > tmp1.csv
sed -e 's/\s\+/,/g' tmp2_fields.csv | sed 's/mRNA/gene/g' | sed 's/\.p1//g' > tmp2.csv

#Retrieve file names
name1=$(dirname $1); name1=$(basename $name1)
name2=$(dirname $2); name2=$(basename $name2)

#Name output directories
outputPath=$(dirname $1); outputPath=$(dirname $outputPath)
outputPath="$outputPath"/Comparisons_GO

#Make output directories
mkdir $outputPath
mkdir "$outputPath"/"$name1"_"$name2"

#Name output files
outputFile="$outputPath"/"$name1"_"$name2"/GO_comparisonSummary.txt
match1File="$outputPath"/"$name1"_"$name2"/"$name1"_GO_matching.csv
unique1File="$outputPath"/"$name1"_"$name2"/"$name1"_GO_unique.csv
match2File="$outputPath"/"$name1"_"$name2"/"$name2"_GO_matching.csv
unique2File="$outputPath"/"$name1"_"$name2"/"$name2"_GO_unique.csv

#Pre-clean up
rm $match1File
rm $unique1File
rm $match2File
rm $unique2File

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp2.csv; then #Matching
		echo "$f1,$f2" >> $match1File
	else #Unique
		echo "$f1,$f2" >> $unique1File
	fi
done < tmp1.csv

#Loop over second set of annotations
while IFS=, read -r f3 f4
do
	#Determine annotation sets
	if grep -q "$f3,$f4" tmp1.csv; then #Matching
		echo "$f3,$f4" >> $match2File
	else #Unique
		echo "$f3,$f4" >> $unique2File
	fi
done < tmp2.csv

#Determine the number of annotations in each subset
total1Num=$(wc -l tmp1.csv | cut -d " " -f 1)
total2Num=$(wc -l tmp2.csv | cut -d " " -f 1)
match1Num=$(wc -l $match1File | cut -d " " -f 1)
match2Num=$(wc -l $match2File | cut -d " " -f 1)
unique1Num=$(wc -l $unique1File | cut -d " " -f 1)
unique2Num=$(wc -l $unique2File | cut -d " " -f 1)

#Output the number of annotations
echo "File Total Unique Matches" > $outputFile
echo "$name1 $total1Num $unique1Num $match1Num" >> $outputFile
echo "$name2 $total2Num $unique2Num $match2Num" >> $outputFile

#Clean up
rm tmp*.csv

echo "GO comparisons complete!"