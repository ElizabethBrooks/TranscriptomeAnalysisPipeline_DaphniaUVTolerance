#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations_GO_trinotate.sh annotationFile1 annotationFile2 annotationFile3
#Usage ex: bash compareAnnotations_GO_trinotate.sh PA42_proteins/user_GO.txt PA42_cds/user_GO.txt PA42_transcripts/user_GO.txt

echo "Beginning GO comparisons..."

#Change commas to semicolons, and white spaces to commas,
# then remove lines without annotations
sed 's/,/;/g' $1 | sed -e 's/\s\+/,/g' | awk -F , 'NF == 2' > tmp1.csv
sed 's/,/;/g' $2 | sed -e 's/\s\+/,/g' | awk -F , 'NF == 2' > tmp2.csv

#Retrieve file names
name1=$(dirname $1); name1=$(basename $name1)
name2=$(dirname $2); name2=$(basename $name2)

#Name output directories
outputPath=$(dirname $1); outputPath=$(dirname $outputPath)
outputPath="$outputPath"/Comparisons_GO_Trinotate

#Pre-clean up
rm -r $outputPath

#Make output directory
mkdir $outputPath

#Name output files
outputFile="$outputPath"/GO_comparisonSummary.txt
file1="$outputPath"/"$name1"_GO.csv
file2="$outputPath"/"$name2"_GO.csv
file1File2="$outputPath"/"$name1"_"$name2"_GO.csv

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp2.csv; then #Intersection 1,2
		echo "$f1,$f2" >> $file1File2
	else #Subset of 1
		echo "$f1,$f2" >> $file1
	fi
done < tmp1.csv

#Loop over second set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp1.csv; then #Intersection 1,2
		echo "$f1,$f2" >> tmp_matched.csv
	else #Subset of 2
		echo "$f1,$f2" >> $file2
	fi
done < tmp2.csv

#Determine the number of annotations in each subset
total1Num=$(wc -l tmp1.csv | cut -d " " -f 1)
total2Num=$(wc -l tmp2.csv | cut -d " " -f 1)
match1Num=$(wc -l $file1 | cut -d " " -f 1)
match2Num=$(wc -l $file2 | cut -d " " -f 1)
match12Num=$(wc -l $file1File2 | cut -d " " -f 1)

#Output the number of annotations
echo "---- Inputs ----" > $outputFile
echo "File Total" >> $outputFile
echo "$name1 $total1Num" >> $outputFile
echo "$name2 $total2Num" >> $outputFile
echo "---- Venn Diagram Sets ----" >> $outputFile
echo "Set Total" >> $outputFile
echo "$name1 $match1Num" >> $outputFile
echo "$name2 $match2Num" >> $outputFile
echo $name1"_"$name2" $match12Num" >> $outputFile

#Clean up
rm tmp*.csv

echo "GO comparisons complete!"