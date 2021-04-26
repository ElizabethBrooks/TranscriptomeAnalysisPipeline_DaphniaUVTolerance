#Retrieve the mapq values for each bam file
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do echo "Processing file $f"; file=$(echo $f | sed "s/accepted_hits.bam/mapq.txt/g"); samtools view $f | cut -f5 > $file; done

#Filter to reads with mapq score of 60 or greater
samtools view -bq 60 /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/accepted_hits.bam  > /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/filtered.bam 

bcftools mpileup -O b -o /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/raw.bcf \
-f /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta \
/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_L2_Pool_1_PA_UV/filtered.bam 

#Keep only unique read alignments
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do echo "Processing file $f"; file=$(echo $f | sed "s/accepted_hits.bam/unique_mapped.bam/g"); samtools view -h -f 0x2 $f | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $file; done

