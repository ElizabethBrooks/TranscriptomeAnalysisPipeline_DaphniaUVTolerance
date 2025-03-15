#Draw UpSet plots of gene sets
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

ddr <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv")
SETDDR <- ddr[,1]

tcast <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/Tcast_Guo_2019_blastp_RBH_PA42_v4.1_geneIDs.csv")
dmel <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/Dmel_Svetec_2016_blastp_RBH_PA42_v4.1_geneIDs.csv")
SETTCAST <- tcast[,1]
SETDMEL <- dmel[,1]

geneCountsInter <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags_filtered.csv")
geneCountsTreat <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv")
geneCountsTol <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags_filtered.csv")

SETInter <- geneCountsInter[,1]
SETTreat <- geneCountsTreat[,1]
SETTol <- geneCountsTol[,1]

lt = list(Interaction=SETInter, Treatment=SETTreat, Tolerance=SETTol, UVT=SETDDR, TCast=SETTCAST, Dmel=SETDMEL)
m = make_comb_mat(lt)

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_combinedDDR_UVR_PA42_v4.1.jpg")
#UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
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

#mUp = m[comb_size(m) >= 20]
#mDown = m[comb_size(m) <= 20]

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_DDRUVR_20Up_PA42_v4.1.jpg")
#UpSet(mUp, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_DDRUVR_20Down_PA42_v4.1.jpg")
#UpSet(mDown, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

lt = list(Interaction=SETInter, Treatment=SETTreat, Tolerance=SETTol, UVT=SETDDR)
m = make_comb_mat(lt)

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_combinedDDR_PA42_v4.1.jpg")
#UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
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

#mUp = m[comb_size(m) >= 20]
#mDown = m[comb_size(m) <= 20]

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_combinedDDR_20Up_PA42_v4.1.jpg")
#UpSet(mUp, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotGLMQLF_combinedDDR_20Down_PA42_v4.1.jpg")
#UpSet(mDown, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

lt = list(Olympics=SETInter, UVT=SETDDR, Tribolium=SETTCAST, Drosophila=SETDMEL)
m = make_comb_mat(lt)

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

#jpeg("/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/plotIntersection_combinedDDR_UVR_PA42_v4.1.jpg")
#UpSet(m, comb_col="#0000FF", bg_col="#F0F0FF", bg_pt_col="#CCCCFF")
#dev.off()

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
