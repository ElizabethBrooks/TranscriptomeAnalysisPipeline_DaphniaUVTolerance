#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations.sh annotationFile1 annotationFile2
#Usage ex: bash compareAnnotations.sh PA42_proteins/user_ko.txt PA42_cds/user_ko.txt

echo "Beginning comparisons!"

#Retrieve file names
name1=$(dirname $1); name1=$(basename $name1)
name2=$(dirname $2); name2=$(basename $name2)

#Name output directory
# and make, if necessary
outputPath=$(dirname $1); outputPath=$(dirname $outputPath)
outputPath="$outputPath"/Comparisons
mkdir $outputPath

#Name output summary file
outputFile="$outputPath"/"$name1"_"$name2"_comparisonSummary.txt

#Change white spaces to commas
# and change mRNA tags to gene
# and remove lines without annotations
sed -e 's/\s\+/,/g' $1 | sed 's/mRNA/gene/g' | sed 's/\.p1//g' | awk -F , 'NF == 2' > tmp1.csv
sed -e 's/\s\+/,/g' $2 | sed 's/mRNA/gene/g' | sed 's/\.p1//g' | awk -F , 'NF == 2' > tmp2.csv

#Clean up
rm "$outputPath"/"$name1"_"$name2"_matching.csv
rm "$outputPath"/"$name1"_unique.csv
rm "$outputPath"/"$name2"_"$name1"_matching.csv
rm "$outputPath"/"$name2"_unique.csv

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp2.csv; then #Matching
		echo "$f1,$f2" >> "$outputPath"/"$name1"_"$name2"_matching.csv
	else #Unique
		echo "$f1,$f2" >> "$outputPath"/"$name1"_unique.csv
	fi
done < tmp1.csv

#Loop over second set of annotations
while IFS=, read -r f3 f4
do
	#Determine annotation sets
	if grep -q "$f3,$f4" tmp1.csv; then #Matching
		echo "$f3,$f4" >> "$outputPath"/"$name2"_"$name1"_matching.csv
	else #Unique
		echo "$f3,$f4" >> "$outputPath"/"$name2"_unique.csv
	fi
done < tmp2.csv

#Determine the number of annotations in each subset
total1Num=$(wc -l tmp1.csv | cut -d " " -f 1)
total2Num=$(wc -l tmp2.csv | cut -d " " -f 1)
matchNum=$(wc -l "$outputPath"/"$name1"_"$name2"_matching.csv | cut -d " " -f 1)
unique1Num=$(wc -l "$outputPath"/"$name1"_unique.csv | cut -d " " -f 1)
unique2Num=$(wc -l "$outputPath"/"$name2"_unique.csv | cut -d " " -f 1)

#Output the number of annotations
echo "---- Inputs ----" > $outputFile
echo "File Total" >> $outputFile
echo "$name1 $total1Num" >> $outputFile
echo "$name2 $total2Num" >> $outputFile
echo "---- Results ---" >> $outputFile
echo "$name1 Unique | Matches | $name2 Unique" >> $outputFile
echo "$unique1Num | $matchNum | $unique2Num" >> $outputFile

#Clean up
rm tmp*.csv

echo "Comparisons complete!"