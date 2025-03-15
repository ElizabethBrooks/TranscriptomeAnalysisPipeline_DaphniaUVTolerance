#Script to retrieve Uniprot IDs from a Trinotate report

cat /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts.csv | cut -f7 | cut -d "^" -f1 > /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/trinotate_annotation_report_PA42_v4.1_transcripts_uniprot.csv
