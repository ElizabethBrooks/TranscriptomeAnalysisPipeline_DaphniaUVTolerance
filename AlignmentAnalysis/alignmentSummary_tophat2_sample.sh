#!/bin/bash
#Bash script to retrieve mapping stats
#Usage: bash alignmentSummary_tophat2_sample.sh alignedSampleFolder analysisMethod runNum
#Usage Ex: bash alignmentSummary_tophat2_sample.sh aligned_tophat2_run2/140327_I481_FCC3P1PACXX_L4_Pool_3_Y05_UV tophat2 2

#Retrieve alignment analysis outputs absolute path
outputsPath=$(grep "alignmentAnalysis:" ../InputData/outputPaths.txt | tr -d " " | sed "s/alignmentAnalysis://g")
#Set outputs directory
outDir="$outputsPath"/AlignmentsAnalyzed
#Prepare input and output file names
runNum=run"$3"
inputStats="$1"align_summary.txt
outputStats="$outDir"/alignmentSummarized_"$2"
#Retrieve sample name
sampleName=$(basename "$1")
#Store sample name in tmp txt file for pasting
echo "$sampleName" > "$outputStats"_sample_"$runNum".txt
#Retrieve mapped left read percents
grep "Mapped" "$inputStats" > "$outputStats"_mapped_"$runNum".txt
head -1 "$outputStats"_mapped_"$runNum".txt > "$outputStats"_mappedLeft_"$runNum".txt
cat "$outputStats"_mappedLeft_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_mappedLeft_"$runNum".txt
sed -i "s/(//g" "$outputStats"_mappedLeft_"$runNum".txt
#Retrieve mapped right read percents
tail -1 "$outputStats"_mapped_"$runNum".txt > "$outputStats"_mappedRight_"$runNum".txt
cat "$outputStats"_mappedRight_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_mappedRight_"$runNum".txt
sed -i "s/(//g" "$outputStats"_mappedRight_"$runNum".txt
#Retrieve overall percents
grep "overall" "$inputStats" > "$outputStats"_overall_"$runNum".txt
cat "$outputStats"_overall_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_overall_"$runNum".txt
#Retrieve concordant percent
grep "concordant" "$inputStats" > "$outputStats"_concordant_"$runNum".txt
cat "$outputStats"_concordant_"$runNum".txt | tr " " "\n" > tmp_"$runNum".txt
grep "%" tmp_"$runNum".txt > "$outputStats"_concordant_"$runNum".txt
#Combine all stats by sample
paste -d "," "$outputStats"_sample_"$runNum".txt "$outputStats"_mappedLeft_"$runNum".txt "$outputStats"_mappedRight_"$runNum".txt "$outputStats"_overall_"$runNum".txt "$outputStats"_concordant_"$runNum".txt > "$outputStats"_combined_"$runNum".csv
#Clean up
rm "tmp_"$runNum".txt"
rm "$outputStats"_sample_"$runNum".txt
rm "$outputStats"_mapped_"$runNum".txt
rm "$outputStats"_mappedLeft_"$runNum".txt
rm "$outputStats"_mappedRight_"$runNum".txt
rm "$outputStats"_overall_"$runNum".txt
rm "$outputStats"_concordant_"$runNum".txt