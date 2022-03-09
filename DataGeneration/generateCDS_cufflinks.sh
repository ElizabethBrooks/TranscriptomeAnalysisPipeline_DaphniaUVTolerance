#!/bin/bash
#Script to generate the longest CDS from the  gff and reference fasta
#Usage: bash generateCDS_cufflinks.sh inputFeatures inputReference genomeTag
#Usage Ex: bash generateCDS_cufflinks.sh /Users/bamflappy/PfrenderLab/NCBI_dataset_Daphnia_pulex_Mar2022/data/KAP4_GCF_021134715.1/genomic.gff /Users/bamflappy/PfrenderLab/NCBI_dataset_Daphnia_pulex_Mar2022/data/KAP4_GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna KAP4
#Usage Ex: bash generateCDS_cufflinks.sh /Users/bamflappy/PfrenderLab/PA42_4.1_resources/PA42.4.1.gff /Users/bamflappy/PfrenderLab/PA42_4.1_resources/GCA_900092285.2_PA42_4.1_genomic.fasta PA42_4.1
#Usage Ex: bash generateCDS_cufflinks.sh /Users/bamflappy/PfrenderLab/PA42_4.2_resources/Genome_assembly_PA42.4.2/PA42.4.2_18449_gene_annotation_file.gff /Users/bamflappy/PfrenderLab/PA42_4.2_resources/Genome_assembly_PA42.4.2/PA42.4.2_genome_assembly.fasta PA42_4.2

#Retrieve input files
inputFeat=$1
inputRef=$2

#Set output file names
outDir=$(dirname $inputRef)
outCDS=$outDir"/"$3"_cds.fa"

#Retrieve CDS
echo "Retrieving CDS..."
#Usage: gffread <input_gff> [-g <genomic_seqs_fasta> | <dir>][-s <seq_info.fsize>] [-o <outfile.gff>] [-t <tname>] [-r #[[<strand>]<chr>:]<start>..<end> [-R]] [-CTVN‐ JMKQAFGUBHZWTOLE] [-w <exons.fa>] [-x <cds.fa>] [-y <tr_cds.fa>] [-i <maxintron>]
gffread "$inputFeat" -g "$inputRef" -x "$outCDS" -W -F