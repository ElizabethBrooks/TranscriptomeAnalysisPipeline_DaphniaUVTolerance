#!/bin/bash
#Script to compare annotation files
#Usage: bash compareAnnotations_KO_ghostkoala.sh annotationFile1 annotationFile2 annotationFile3
#Usage ex: bash compareAnnotations_KO_ghostkoala.sh PA42_proteins/user_ko.txt PA42_cds/user_ko.txt PA42_transcripts/user_ko.txt

echo "Beginning KO comparisons..."

#Change white spaces to commas
# and change mRNA tags to gene
# and remove lines without annotations
sed -e 's/\s\+/,/g' $1 | sed 's/mRNA/gene/g' | sed 's/\.p.//g' | awk -F , 'NF == 2' > tmp1.csv
sed -e 's/\s\+/,/g' $2 | sed 's/mRNA/gene/g' | sed 's/\.p.//g' | awk -F , 'NF == 2' > tmp2.csv
sed -e 's/\s\+/,/g' $3 | sed 's/mRNA/gene/g' | sed 's/\.p.//g' | awk -F , 'NF == 2' > tmp3.csv

#Retrieve file names
name1=$(dirname $1); name1=$(basename $name1)
name2=$(dirname $2); name2=$(basename $name2)
name3=$(dirname $3); name3=$(basename $name3)

#Name output directories
outputPath=$(dirname $1); outputPath=$(dirname $outputPath)
outputPath="$outputPath"/Comparisons_KO/"$name1"_"$name2"_"$name3"

#Pre-clean up
rm -r $outputPath

#Make output directory
mkdir $outputPath

#Name output files
outputFile="$outputPath"/KO_comparisonSummary.txt
file1="$outputPath"/"$name1"_KO.csv
file2="$outputPath"/"$name2"_KO.csv
file3="$outputPath"/"$name3"_KO.csv
file1File2="$outputPath"/"$name1"_"$name2"_KO.csv
file1File3="$outputPath"/"$name1"_"$name3"_KO.csv
file2File3="$outputPath"/"$name2"_"$name3"_KO.csv
file1File2File3="$outputPath"/"$name1"_"$name2"_"$name3"_KO.csv

#Loop over first set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp2.csv; then
		if grep -q "$f1,$f2" tmp3.csv; then #Intersection 1,2,3
			echo "$f1,$f2" >> $file1File2File3
		else #Intersection 1,2
			echo "$f1,$f2" >> $file1File2
		fi
	elif grep -q "$f1,$f2" tmp3.csv; then  #Intersection 1,3
		echo "$f1,$f2" >> $file1File3
	else #Subset of 1
		echo "$f1,$f2" >> $file1
	fi
done < tmp1.csv

#Loop over second set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp1.csv; then 
		if grep -q "$f1,$f2" tmp3.csv; then #Intersection 1,2,3
			echo "$f1,$f2" >> tmp2_matches.csv
		else #Intersection 1,2
			echo "$f1,$f2" >> tmp2_matches.csv
		fi
	elif grep -q "$f1,$f2" tmp3.csv; then  #Intersection 2,3
		echo "$f1,$f2" >> $file2File3
	else #Subset of 2
		echo "$f1,$f2" >> $file2
	fi
done < tmp2.csv

#Loop over third set of annotations
while IFS=, read -r f1 f2
do
	#Determine annotation sets
	if grep -q "$f1,$f2" tmp1.csv; then 
		if grep -q "$f1,$f2" tmp2.csv; then #Intersection 1,2,3
			echo "$f1,$f2" >> tmp3_matches.csv
		else #Intersection 1,3
			echo "$f1,$f2" >> tmp3_matches.csv
		fi
	elif grep -q "$f1,$f2" tmp2.csv; then  #Intersection 2,3
		echo "$f1,$f2" >> tmp3_matches.csv
	else #Subset of 3
		echo "$f1,$f2" >> $file3
	fi
done < tmp3.csv

#Determine the number of annotations in each subset
total1Num=$(wc -l tmp1.csv | cut -d " " -f 1)
total2Num=$(wc -l tmp2.csv | cut -d " " -f 1)
total3Num=$(wc -l tmp3.csv | cut -d " " -f 1)
match1Num=$(wc -l $file1 | cut -d " " -f 1)
match2Num=$(wc -l $file2 | cut -d " " -f 1)
match3Num=$(wc -l $file3 | cut -d " " -f 1)
match12Num=$(wc -l $file1File2 | cut -d " " -f 1)
match23Num=$(wc -l $file2File3 | cut -d " " -f 1)
match13Num=$(wc -l $file1File3 | cut -d " " -f 1)
match123Num=$(wc -l $file1File2File3 | cut -d " " -f 1)

#Output the number of annotations
echo "---- Inputs ----" > $outputFile
echo "File Total" >> $outputFile
echo "$name1 $total1Num" >> $outputFile
echo "$name2 $total2Num" >> $outputFile
echo "$name3 $total3Num" >> $outputFile
echo "---- Venn Diagram Sets ----" >> $outputFile
echo "Set Total" >> $outputFile
echo "$name1 $match1Num" >> $outputFile
echo "$name2 $match2Num" >> $outputFile
echo "$name3 $match3Num" >> $outputFile
echo $name1"_"$name2" $match12Num" >> $outputFile
echo $name1"_"$name3" $match13Num" >> $outputFile
echo $name2"_"$name3" $match23Num" >> $outputFile
echo $name1"_"$name2"_"$name3" $match123Num" >> $outputFile

#Clean up
rm tmp*.csv

echo "KO comparisons complete!"