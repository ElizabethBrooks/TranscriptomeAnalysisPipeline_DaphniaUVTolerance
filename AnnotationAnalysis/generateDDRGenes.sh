dir="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1"

grep "GO:1990391" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 > "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0031570" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0009411" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0007093" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0006974" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0003697" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
grep "GO:0003684" "$dir"/gene2GO_PA42_v4.1_transcripts.map | cut -d$'\t' -f 1 >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
echo "dp_gene2015" >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
echo "dp_gene2349" >> "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv
cat "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv | awk '!seen[$0] {print} {++seen[$0]}' > "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs_uniq.csv
rm "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs.csv

cat "$dir"/DDR_Dmel_Svetec_2016_blastp_RBH_PA42_v4.1_geneIDs.csv > "$dir"/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs.csv
cat "$dir"/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs_uniq.csv >> "$dir"/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs.csv
cat "$dir"/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs.csv | awk '!seen[$0] {print} {++seen[$0]}' > "$dir"/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv
rm "$dir"/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs.csv