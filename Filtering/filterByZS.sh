#Keep only unique read alignments using ZS tag
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do echo "Processing file $f"; path=$(dirname $f); samtools view -h -f 0x2 $f | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > "$path"/filteredZS.bam; done
