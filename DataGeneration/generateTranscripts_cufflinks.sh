#!/bin/bash
#Script to generate a multi FASTA file for all transcripts in a GFF file
#Usage: bash generateTranscripts_cufflinks.sh genomeFeat genomeRef
#Usage Ex: bash generateTranscripts_cufflinks.sh /Users/bamflappy/PfrenderLab/NCBI_dataset_Daphnia_pulex_Mar2022/data/KAP4_GCF_021134715.1/genomic.gff /Users/bamflappy/PfrenderLab/NCBI_dataset_Daphnia_pulex_Mar2022/data/KAP4_GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna KAP4
#Usage Ex: bash generateTranscripts_cufflinks.sh /Users/bamflappy/PfrenderLab/PA42_4.1_resources/PA42.4.1.gff /Users/bamflappy/PfrenderLab/PA42_4.1_resources/GCA_900092285.2_PA42_4.1_genomic.fasta PA42_4.1
#Usage Ex: bash generateTranscripts_cufflinks.sh /Users/bamflappy/PfrenderLab/PA42_4.2_resources/Genome_assembly_PA42.4.2/PA42.4.2_18449_gene_annotation_file.gff /Users/bamflappy/PfrenderLab/PA42_4.2_resources/Genome_assembly_PA42.4.2/PA42.4.2_genome_assembly.fasta PA42_4.2

#Retrieve genome reference and features absolute paths
inputFeatFile=$1
inputRef=$2

#Name output file of inputs
outputsPath=$(dirname $inputRef)
inputOutFile=$outputsPath"/"$3"_generateTranscipts_summary.txt"
errorOut=$outputsPath"/"$3"_generateTransciptsError_summary.txt"

#Generate a fasta index file using samtools
samtools faidx $inputRef

#Generate a FASTA file with the DNA sequences for all transcripts in the GFF file
gffread -w $outputsPath"/"$3"_transcripts_cufflinks.fa" -g $inputRef $inputFeatFile 2> "$errorOut"
echo "Generate fasta file with the DNA sequences for all transcripts in the GFF file..." > "$inputOutFile"
echo "gffread -w $outputsPath/$3_transcripts_cufflinks.fa -g $inputRef $inputFeatFile" >> "$inputOutFile"
