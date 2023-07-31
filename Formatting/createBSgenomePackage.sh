

# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

cat /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna | grep ">" | cut -d " " -f 1 | sed "s/>//g" > tmp.txt

while read -r line; do cat /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/GCF_021134715.1_ASM2113471v1_genomic.fna | tr -s '\n' ';' | tr -s '>' '\n' | grep $line | sed "s/^/>/g" | tr -s ';' '\n' > /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs/"$line".fa ; done < tmp.txt

#library(BSgenome)
#setwd("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs")
#forgeBSgenomeDataPkg("/Users/bamflappy/Repos/TranscriptomeAnalysisPipeline_DaphniaUVTolerance/InputData/BSgenome.Dpulex.KAP4-seed", replace=TRUE)

R CMD build /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs/BSgenome.Dpulex.KAP4

R CMD check /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs/BSgenome.Dpulex.KAP4_1.0.0.tar.gz

R CMD INSTALL /Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/ncbi_dataset/data/GCF_021134715.1/seqs/BSgenome.Dpulex.KAP4_1.0.0.tar.gz
