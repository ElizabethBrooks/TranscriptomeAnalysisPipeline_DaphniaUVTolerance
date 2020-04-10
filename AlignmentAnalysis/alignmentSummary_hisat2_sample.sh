#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_hisat2_sample.sh alignedSampleFolder alignmentMethod runNum
#Usage Ex: bash alignmentSummary_hisat2_sample.sh aligned_hisat2_run1/140327_I481_FCC3P1PACXX_L4_Pool_3_Y05_UV hisat2 run1
#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory
outDir="$outputsPath"/AlignmentsAnalyzed
#Prepare input and output file names
runNum="$3"
inputStats="$1"/alignedSummary.txt
outputStats="$outDir"/alignmentSummarized_"$2"
#Retrieve sample name
sampleName=$(basename "$1")
#Store sample name in tmp txt file for pasting
echo "$sampleName" > "$outputStats"_sample_"$runNum".txt
#Retrieve concordant alignment percents
grep "aligned concordantly exactly" "$inputStats" > "$outputStats"_concordant_"$runNum".txt
cat "$outputStats"_concordant_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_concordant_"$runNum".txt
sed -i "s/(//g" "$outputStats"_concordant_"$runNum".txt
sed -i "s/)//g" "$outputStats"_concordant_"$runNum".txt
#Retrieve overall alignment percent
grep "overall alignment rate" "$inputStats" > "$outputStats"_overall_"$runNum".txt
cat "$outputStats"_overall_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_overall_"$runNum".txt
#Combine all stats by sample
paste -d "," "$outputStats"_sample_"$runNum".txt "$outputStats"_overall_"$runNum".txt "$outputStats"_concordant_"$runNum".txt > "$outputStats"_combined_"$runNum".csv
#Clean up
rm "tmp_"$runNum".txt"
rm "$outputStats"_sample_"$runNum".txt
rm "$outputStats"_concordant_"$runNum".txt
rm "$outputStats"_overall_"$runNum".txt