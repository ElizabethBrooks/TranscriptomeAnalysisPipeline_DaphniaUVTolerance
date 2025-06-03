#Script to map gene IDs to Uniprot IDs
#Usage: bash generateGeneUniprotMap.sh

#Retrieve gene IDs
cat /Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts.csv | cut -d$'\t' -f1 > tmp1.txt

#Retrieve Uniprot IDs
cat /Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts.csv | cut -d$'\t' -f7 | cut -d "^" -f1 | cut -d "_" -f1 > tmp2.txt

#Create a mapping file of gene and Uniprot IDs
paste -d "," tmp1.txt tmp2.txt > tmp.csv

#Squeeze any sequence of newline characters to one
cat tmp.csv | tr -s '\n' > /Users/bamflappy/PfrenderLab/PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_geneUniprotMap.csv

#Clean up
rm tmp*