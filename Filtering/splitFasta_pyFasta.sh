#!/bin/bash
#Script to split a fasta file into relatively even pieces,
# while maintaining whole sequences
#Usage: bash splitFasta_pyFasta.sh numPieces
#Alternate usage: bash splitFasta_pyFasta.sh 4

#Split fasta into the user input number of pieces,
# keeping whole seqs
pyfasta split -n $1 Daphnia_pulex_PA42_v3.0.fasta

#Check number of sequences
totalSeqs=0
echo "Number of sequences:"
for f in Daphnia_pulex_PA42_v3.0.*.fasta; do
	seqs=$(grep ">" $f | wc -l); echo "$seqs $f"
	totalSeqs=$(($totalSeqs + $seqs))
done
echo "$totalSeqs total"

#Check number of lines
echo "Number of lines:"
wc -l Daphnia_pulex_PA42_v3.0.*.fasta

#Check file sizes
echo "File sizes (bytes):"
wc -c Daphnia_pulex_PA42_v3.0.*.fasta