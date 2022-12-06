#!/bin/bash
#Script to test the PA42 reference genome fasta formatting
#Usage: bash splitFasta_test.sh userName localDir
#Usage ex: bash splitFasta_test.sh ebrooks5 ~/Documents/

#Copy PA42 reference genome fasta
scp "$1"@crcfe01.crc.nd.edu:/afs/crc.nd.edu/group/pfrenderlab/mendel/ebrooks/rnaseq/Daphnia_pulex_PA42_v3.0.fasta "$2"

#Remove all newlines from lines that do not start with the fasta header character
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' Daphnia_pulex_PA42_v3.0.fasta > Daphnia_pulex_PA42_v3.0_cleaned.fasta

#Convert a multiline fasta to singleline
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Daphnia_pulex_PA42_v3.0_cleaned.fasta > Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta

#Split a fasta file
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("myseq_singleline%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta 

#Check number of sequences
echo "Number of sequences:"
seqs=$(grep ">" Daphnia_pulex_PA42_v3.0.fasta | wc -l); echo "$seqs Daphnia_pulex_PA42_v3.0.fasta"
seqs=$(grep ">" Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta | wc -l); echo "$seqs Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta"
seqs=$(grep ">" myseq_singleline0.fa | wc -l); echo "$seqs myseq_singleline0.fa"
seqs=$(grep ">" myseq_singleline1000.fa | wc -l); echo "$seqs myseq_singleline1000.fa"

#Check number of lines
echo "Number of lines:"
wc -l Daphnia_pulex_PA42_v3.0.fasta
wc -l Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta
wc -l myseq_singleline0.fa
wc -l myseq_singleline1000.fa

#Check file sizes
echo "File sizes (bytes):"
wc -c Daphnia_pulex_PA42_v3.0.fasta
wc -c Daphnia_pulex_PA42_v3.0_cleaned_singleline.fasta
wc -c myseq_singleline0.fa
wc -c myseq_singleline1000.fa

#Remove duplicate seqs in a fasta file
#Note that seqkit will need to be installed first
#cat Daphnia_pulex_PA42_v3.0.fasta | seqkit rmdup -s -i -o clean.fasta -d duplicated.fasta -D duplicated.detail.txt