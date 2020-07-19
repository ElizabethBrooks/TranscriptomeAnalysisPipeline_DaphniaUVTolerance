#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations.sh annotationType annotationFile1 annotationFile2
#Usage ex: bash compareAnnotations.sh KO PA42_proteins/user_ko.txt PA42_cds/user_ko.txt
#Usage ex: bash compareAnnotations.sh GO PA42_proteins/user_ko.txt PA42_cds/user_ko.txt

echo "Beginning comparisons..."

#Retrieve file names
name1=$(dirname $2); name1=$(basename $name1)
name2=$(dirname $3); name2=$(basename $name2)

#Change white spaces to commas
# and change mRNA tags to gene
# and remove lines without annotations
sed -e 's/\s\+/,/g' $2 | sed 's/mRNA/gene/g' | sed 's/\.p1//g' | awk -F , 'NF == 2' > tmp1.csv
sed -e 's/\s\+/,/g' $3 | sed 's/mRNA/gene/g' | sed 's/\.p1//g' | awk -F , 'NF == 2' > tmp2.csv

#Name output directory
# and make, if necessary
outputPath=$(dirname $2); outputPath=$(dirname $outputPath)
outputPath="$outputPath"/Comparisons_"$1"
mkdir $outputPath

#Name output files
outputFile="$outputPath"/"$name1"_"$name2"_"$1"_comparisonSummary.txt
match1File="$outputPath"/"$name1"_"$name2"_"$1"_matching.csv
unique1File="$outputPath"/"$name1"_"$1"_unique.csv
match2File="$outputPath"/"$name2"_"$name1"_"$1"_matching.csv
unique2File="$outputPath"/"$name2"_"$1"_unique.csv

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
matchNum=$(wc -l $match1File | cut -d " " -f 1)
unique1Num=$(wc -l $unique1File | cut -d " " -f 1)
unique2Num=$(wc -l $unique2File | cut -d " " -f 1)

#Output the number of annotations
echo "---- $1 Inputs ----" > $outputFile
echo "File Total" >> $outputFile
echo "$name1 $total1Num" >> $outputFile
echo "$name2 $total2Num" >> $outputFile
echo "---- $1 Comparison ---" >> $outputFile
echo "$name1 | Matches | $name2" >> $outputFile
echo "$unique1Num | $matchNum | $unique2Num" >> $outputFile

#Clean up
rm tmp*.csv

echo "Comparisons complete!"