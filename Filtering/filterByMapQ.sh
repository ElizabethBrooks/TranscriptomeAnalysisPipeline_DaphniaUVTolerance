#Keep only unique read alignments usingthe mapq score of 60 
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do echo "Processing file $f";  path=$(dirname $f); samtools view -bq 60 $f > "$path"/filteredMapQ.bam; done
