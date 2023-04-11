#!/bin/bash
#$ -M ebrooks5@nd.edu
#$ -m abe
#$ -r n
#$ -N translateCodons_jobOutput

# script to generate codon alignments for each gene in the reference set of peptide sequences
# usage: bash translateCodons_pal2nal.sh
# usage ex: bash translateCodons_pal2nal.sh


#Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]
echo "Generating codon alignment for $gTag..."
#"$softwarePath"/pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  "$outPath"/"$gTag".codon
pal2nal.pl "$inAln" "$tmpRefNuc" "$tmpConNuc" -output paml -nogap  >  "$outPath"/"$gTag".codon
