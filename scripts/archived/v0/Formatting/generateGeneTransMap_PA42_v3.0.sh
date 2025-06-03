#!/bin/bash

grep ">" /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/Frozen3.0Genome_Annotaton/PA42.3.0.protein_new.fasta | sed "s/>//g" > /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col1.txt
tr -d "\n\r" < /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col1.txt | sed "s/gene/\ngene/g" > /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col1_cleaned.txt

grep ">" /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/Frozen3.0Genome_Annotaton/PA42.3.0.cds_new.fasta | sed "s/>//g" > /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col2.txt
tr -d "\n\r" < /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col2.txt | sed "s/mRNA/\nmRNA/g" > /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col2_cleaned.txt

paste /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col1_cleaned.txt /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/col2_cleaned.txt > /home/mae/Documents/RNASeq_Workshop_ND/PA42_annotations/PA42.3.0.fasta.gene_trans_map