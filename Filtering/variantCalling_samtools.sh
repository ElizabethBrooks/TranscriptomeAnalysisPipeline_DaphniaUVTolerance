#!/bin/bash
#Script to perform variant calling

#Genome reference file
ref = "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta"

#Loop over MapQ filtered bam files
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/filteredMapQ.bam; do 
	echo "Processing file $f"
	path=$(dirname $f)
	file=$(basename $f | sed "s/.bam//g")
	#Calculate the read coverage of positions in the genome
	bcftools mpileup -Ob -o "$path"/"$file"_raw.bcf -f "/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta" "$f" 
	#Detect the single nucleotide polymorphisms 
	bcftools call -mv -Ob -o "$path"/"$file"_calls.vcf.gz "$path"/"$file"_raw.bcf 
	#Index vcf file
	bcftools index "$path"/"$file"_calls.vcf.gz
	#Normalize indels
	#bcftools norm -f "$ref" "$path"/"$file"_calls.vcf.gz -Ob -o "$path"/"$file"_calls.norm.bcf
	#Filter adjacent indels within 5bp
	#bcftools filter --IndelGap 5 "$path"/"$file"_calls.norm.bcf -Ob -o "$path"/"$file"_calls.norm.flt-indels.bcf
	#Exclude sites where FILTER is not set
	#bcftools query -e'FILTER="."' -f'%CHROM %POS %FILTER\n' "$path"/"$file"_calls.norm.flt-indels.bcf > "$path"/"$file"_filtered_excluded.bcf
done