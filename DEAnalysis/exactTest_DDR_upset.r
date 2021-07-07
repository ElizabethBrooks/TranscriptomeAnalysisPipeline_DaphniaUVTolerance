#Draw UpSet plots of gene sets

#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

ddr <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv")
SETDDR <- ddr[,1]

geneCountsE05 <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestE05Analysis_FDR0.10/exactTest_topTags_filtered.csv")
geneCountsY05 <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestY05Analysis_FDR0.10/exactTest_topTags_filtered.csv")
geneCountsY023 <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestY023_5Analysis_FDR0.10/exactTest_topTags_filtered.csv")
geneCountsR2 <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestR2Analysis_FDR0.10/exactTest_topTags_filtered.csv")
geneCountsPA <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestPAAnalysis_FDR0.10/exactTest_topTags_filtered.csv")
geneCountsSierra <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/DEAnalysis_run1/exactTestSierraAnalysis_FDR0.10/exactTest_topTags_filtered.csv")

SETE05 <- geneCountsE05[,1] 
SETY05 <- geneCountsY05[,1]
SETY023 <- geneCountsY023[,1]
SETR2 <- geneCountsR2[,1]
SETPA <- geneCountsPA[,1]
SETSierra <- geneCountsSierra[,1]

SETE05 <- SETE05[SETE05 %in% SETDDR]
SETY05 <- SETY05[SETY05 %in% SETDDR]
SETY023 <- SETY023[SETY023 %in% SETDDR]
SETR2 <- SETR2[SETR2 %in% SETDDR]
SETPA <- SETPA[SETPA %in% SETDDR]
SETSierra <- SETSierra[SETSierra %in% SETDDR]

lt = list(E05=SETE05, Y05=SETY05, Y023=SETY023, R2=SETR2, PA=SETPA, Sierra=SETSierra)
m = make_comb_mat(lt)

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/plotExactTests_combinedDDR_PA42_v4.1.jpg")
#UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

#sizeUp = m[comb_size(m) >= 10]
#sizeDown = m[comb_size(m) <= 10]

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/plotExactTests_combinedDDR_20Up_PA42_v4.1.jpg")
#UpSet(sizeUp, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/plotExactTests_combinedDDR_20Down_PA42_v4.1.jpg")
#UpSet(sizeDown, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

cs = comb_size(m)
ht = UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
plot <- decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 5), 
            just = "bottom",
            default.units = "native")
})

size = m[comb_degree(m) >= 2]
cs = comb_size(size)
ht = UpSet(size, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF", top_annotation = upset_top_annotation(size, ylim = c(0, 1.1*max(cs))))
ht = draw(ht)
co = column_order(ht)
nc = ncol(size)
plot <- decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 5), 
            just = "bottom",
            default.units = "native")
})

#m1 = make_comb_mat(lt, mode = "distinct")
#m2 = make_comb_mat(lt, mode = "intersect")
#m3 = make_comb_mat(lt, mode = "union")
#UpSet(m1, row_title = "distinct mode")
#UpSet(m2, row_title = "intersect mode")
#UpSet(m3, row_title = "union mode")