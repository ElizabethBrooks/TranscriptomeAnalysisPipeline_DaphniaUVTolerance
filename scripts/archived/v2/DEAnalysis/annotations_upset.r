#Draw UpSet plots of gene sets
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

geneList <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/PA42.4.1_geneIDs.csv")
ddrPA <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_PA42_v4.1_transcripts.map_geneIDs_uniq.csv")
ddrDro <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDR_Dmel_Svetec_2016_blastp_RBH_PA42_v4.1_geneIDs.csv")

SETPA <- geneList[,1]
SETDDRPA <- ddrPA[,1]
SETDDRDro <- ddrDro[,1]

lt = list(DaphniaGenes=SETPA, DaphniaUVT=SETDDRPA, DrosophilaUVT=SETDDRDro)
m = make_comb_mat(lt)

cs = comb_size(m)
ht = UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
pdf("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_nums.pdf")
plot <- decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 5), 
            just = "bottom",
            default.units = "native")
})
dev.off()

size = m[comb_degree(m) >= 2]
cs = comb_size(size)
ht = UpSet(size, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", top_annotation = upset_top_annotation(size, ylim = c(0, 1.1*max(cs))))
ht = draw(ht)
co = column_order(ht)
nc = ncol(size)
pdf("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_2UpNums.pdf")
plot <- decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 5), 
            just = "bottom",
            default.units = "native")
})
dev.off()

#m1 = make_comb_mat(lt, mode = "distinct")
#m2 = make_comb_mat(lt, mode = "intersect")
#m3 = make_comb_mat(lt, mode = "union")
#UpSet(m1, row_title = "distinct mode")
#UpSet(m2, row_title = "intersect mode")
#UpSet(m3, row_title = "union mode")
