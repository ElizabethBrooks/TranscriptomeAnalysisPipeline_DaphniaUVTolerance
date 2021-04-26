#Perform variant calling
for f in /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/sortedCoordinate_samtoolsHisat2_run2/140327_I481_FCC3P1PACXX_*/accepted_hits.bam; do 
	echo "Processing file $f"
	path=$(dirname $f)
	#Calculate the read coverage of positions in the genome
	bcftools mpileup -O b -o "$path"/raw.bcf -f /home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/GCA_900092285.2_PA42_4.1_genomic.fasta "$path"/filteredMapQ.bam 
	#Detect the single nucleotide polymorphisms 
	bcftools call --ploidy 2 -m -v -o "$path"/variants.vcf "$path"/raw.bcf 
	#Filter and report the SNP variants in variant calling format
	vcfutils.pl varFilter "$path"/variants.vcf > "$path"/final_variants.vcf
	#Count the number of variants
	grep -v "#" "$path"/final_variants.vcf | wc -l
done