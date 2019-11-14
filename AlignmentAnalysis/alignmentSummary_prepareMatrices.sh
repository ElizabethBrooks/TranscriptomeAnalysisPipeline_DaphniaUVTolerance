#!/bin/bash
#Script to prepare alignment summary matrices
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