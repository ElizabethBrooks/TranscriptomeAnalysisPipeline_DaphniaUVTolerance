#!/bin/bash
#Script to prepare alignment summary matrices

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Move to outputs directory
cd "$outputsPath"/AlignmentsAnalyzed
#Create directory for alignment analysis
#outputAnalysis=AlignmentAnalysis
#mkdir "$outputAnalysis"
#Move to folder with alignment summaries
cd AlignmentStats_Analysis

sed -i "s/%//g" alignmentSummarized_tophat2.csv 

sed -i "s/140327_I481_FCC3P1PACXX_L2_//g" alignmentSummarized_tophat2_run2.csv
sed -i "s/140327_I481_FCC3P1PACXX_L3_//g" alignmentSummarized_tophat2_run2.csv
sed -i "s/140327_I481_FCC3P1PACXX_L4_//g" alignmentSummarized_tophat2_run2.csv

sed -i "s/Pool_1/Pool1/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/Pool_2/Pool2/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/Pool_3/Pool3/g" alignmentSummarized_tophat2_run2.csv

sed -i "s/UV1/UV/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/UV2/UV/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/UV3/UV/g" alignmentSummarized_tophat2_run2.csv

sed -i "s/VIS1/UV/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/VIS2/UV/g" alignmentSummarized_tophat2_run2.csv
sed -i "s/VIS3/UV/g" alignmentSummarized_tophat2_run2.csv

head -1 alignmentSummarized_tophat2_run2.csv > tmp.csv

grep "Pool1.*UV" alignmentSummarized_tophat2_run2.csv | sed "s/UV/UV_Pool1/g" >> tmp.csv
grep "Pool2.*UV" alignmentSummarized_tophat2_run2.csv | sed "s/UV/UV_Pool2/g" >> tmp.csv
grep "Pool3.*UV" alignmentSummarized_tophat2_run2.csv | sed "s/UV/UV_Pool3/g" >> tmp.csv

grep "Pool1.*VIS" alignmentSummarized_tophat2_run2.csv | sed "s/VIS/VIS_Pool1/g" >> tmp.csv
grep "Pool2.*VIS" alignmentSummarized_tophat2_run2.csv | sed "s/VIS/VIS_Pool2/g" >> tmp.csv
grep "Pool3.*VIS" alignmentSummarized_tophat2_run2.csv | sed "s/VIS/VIS_Pool3/g" >> tmp.csv

sed -i "s/Pool1_//g" tmp.csv
sed -i "s/Pool2_//g" tmp.csv
sed -i "s/Pool3_//g" tmp.csv

mv tmp.csv alignmentSummarized_tophat2.csv

head -1 alignmentSummarized_tophat2.csv > tmp.csv
tail -n +2 alignmentSummarized_tophat2.csv | sort -k1 -n -t, >> tmp.csv

mv tmp.csv alignmentSummarized_tophat2.csv

cut -d, -f2-3 --complement alignmentSummarized_tophat2.csv > alignmentSummarized_tophat2_subset.csv 

sed 's/$/,tophat2/' alignmentSummarized_tophat2.csv > alignmentSummarized_tophat2_tagged.csv 
sed -i 's/concordant,tophat2/concordant,method/' alignmentSummarized_tophat2_tagged.csv

cut -d, -f1 --complement alignmentSummarized_tophat2_subset.csv > tmp.csv
sed "s/overall/overallTophat2/g" tmp.csv > alignmentSummarized_tophat2_subset_trimmed.csv
sed -i "s/concordant/concordantTophat2/g" alignmentSummarized_tophat2_subset_trimmed.csv

cut -d, -f1 alignmentSummarized_legacy_subset_trimmed.csv > tmpL1.csv
cut -d, -f1 alignmentSummarized_tophat2_subset_trimmed.csv > tmpT1.csv
cut -d, -f1 alignmentSummarized_hisat2_trimmed.csv > tmpH1.csv
cut -d, -f2 alignmentSummarized_legacy_subset_trimmed.csv > tmpL2.csv
cut -d, -f2 alignmentSummarized_tophat2_subset_trimmed.csv > tmpT2.csv
cut -d, -f2 alignmentSummarized_hisat2_trimmed.csv > tmpH2.csv
paste -d "," tmpL1.csv tmpT1.csv tmpH1.csv tmpL2.csv tmpT2.csv tmpH2.csv > alignmentSummarized_legacyTophat2Hisat2_differences_merged.csv
rm tmp*.csv

bash generateDifferences.sh alignmentSummarized_legacy_subset_trimmed.csv alignmentSummarized_tophat2_subset_trimmed.csv

paste -d "," alignmentSummarized_legacyTophat2_differences_overall.csv alignmentSummarized_legacyTophat2_differences_concordant.csv > alignmentSummarized_legacyTophat2_differences.csv
sed 's/$/,tophat2/' alignmentSummarized_legacyTophat2_differences.csv > tmp.csv 
sed -i 's/concordantDifferences,tophat2/concordantDifferences,method/' tmp.csv

cut -d, -f1 alignmentSummarized_tophat2.csv > samples.csv
paste -d "," samples.csv tmp.csv > alignmentSummarized_legacyTophat2_differences_tagged.csv
rm tmp.csv
rm samples.csv

echo "sampleNum" > lineNums.csv
for i in {1..36}; do echo "$i" >> lineNums.csv; done
paste -d "," lineNums.csv alignmentSummarized_legacyTophat2_differences_tagged.csv >> alignmentSummarized_legacyTophat2_differences_numbered.csv
paste -d "," lineNums.csv alignmentSummarized_legacyHisat2_differences_tagged.csv >> alignmentSummarized_legacyHisat2_differences_numbered.csv
rm lineNums.csv

grep -v "sampleNum" alignmentSummarized_legacyHisat2_differences_numbered.csv > tmp.csv
cat alignmentSummarized_legacyTophat2_differences_numbered.csv > alignmentSummarized_legacyTophat2Hisat2_differences_numbered.csv
cat tmp.csv >> alignmentSummarized_legacyTophat2Hisat2_differences_numbered.csv
rm tmp.csv

echo "line" > lineNums.csv
for i in {1..72}; do echo "$i" >> lineNums.csv; done
paste -d "," lineNums.csv alignmentSummarized_legacyTophat2Hisat2_differences_numbered.csv >> alignmentSummarized_legacyTophat2Hisat2_differences_merged.csv
rm lineNums.csv

sed 's/$/,legacy/' alignmentSummarized_legacy_subset.csv > tmp.csv 
sed 's/concordant,legacy/concordant,method/' tmp.csv > alignmentSummarized_legacy_subset_tagged.csv
rm tmp.csv

echo "sampleNum" > lineNums.csv
for i in {1..36}; do echo "$i" >> lineNums.csv; done
paste -d "," lineNums.csv alignmentSummarized_legacy_subset_tagged.csv >> alignmentSummarized_legacy_subset_numbered.csv
paste -d "," lineNums.csv alignmentSummarized_tophat2_subset_tagged.csv >> alignmentSummarized_tophat2_subset_numbered.csv
paste -d "," lineNums.csv alignmentSummarized_hisat2_tagged.csv >> alignmentSummarized_hisat2_numbered.csv
rm lineNums.csv

grep -v "sampleNum" alignmentSummarized_tophat2_subset_numbered.csv > tmp1.csv
grep -v "sampleNum" alignmentSummarized_hisat2_numbered.csv > tmp2.csv
cat alignmentSummarized_legacy_subset_numbered.csv > alignmentSummarized_legacyTophat2Hisat2_numbered.csv
cat tmp1.csv >> alignmentSummarized_legacyTophat2Hisat2_numbered.csv
cat tmp2.csv >> alignmentSummarized_legacyTophat2Hisat2_numbered.csv
rm tmp*

echo "line" > lineNums.csv
for i in {1..108}; do echo "$i" >> lineNums.csv; done
paste -d "," lineNums.csv alignmentSummarized_legacyTophat2Hisat2_numbered.csv >> alignmentSummarized_legacyTophat2Hisat2_merged.csv
rm lineNums.csv