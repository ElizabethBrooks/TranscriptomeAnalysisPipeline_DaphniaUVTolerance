#Retrieve the gene IDs
grep ">" /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/transcripts_cufflinks.fa | sed "s/>//g" | cut -d" " -f1 | cut -d"-" -f1 > col1.txt

#Retrieve the gtranscript IDs
grep ">" /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/transcripts_cufflinks.fa | sed "s/>//g" | cut -d" " -f1 > col2.txt

#Combine IDs to make the gene to transcript map
paste -d "\t" col1.txt col2.txt > /afs/crc.nd.edu/group/pfrenderlab/bateson/ebrooks/rnaseq/genomicResources_PA42_v4.1/transcripts_cufflinks.fa.gene_trans_map

#Clean up
rm col*.txt