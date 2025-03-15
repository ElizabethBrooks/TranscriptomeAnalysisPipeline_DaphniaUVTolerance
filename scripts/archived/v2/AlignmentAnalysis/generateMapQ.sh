#Retrieve the mapq values for each bam file
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do echo "Processing file $f"; file=$(echo $f | sed "s/accepted_hits.bam/mapq.txt/g"); samtools view $f | cut -f5 > $file; done
