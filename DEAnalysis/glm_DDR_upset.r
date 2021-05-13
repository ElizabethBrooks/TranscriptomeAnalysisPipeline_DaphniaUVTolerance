#Draw UpSet plots of gene sets
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

ddr <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv")
SETDDR <- ddr[,1]

tcast <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/Tcast_Guo_2019_blastp_RBH_PA42_v4.1_geneIDs.csv")
dmel <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/Dmel_Svetec_2016_blastp_RBH_PA42_v4.1_geneIDs.csv")
INTCAST <- tcast[,1]
INDMEL <- dmel[,1]

SETTCAST <- INTCAST[INTCAST %in% SETDDR]
SETDMEL <- INDMEL[INDMEL %in% SETDDR]

geneCountsInter <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_interaction_topTags_filtered.csv")
geneCountsTreat <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_UVvsVIS_topTags_filtered.csv")
geneCountsTol <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/genomicResources_PA42_v4.1/geneCounts_mergedHisat2_PA42_v4.1/glmQLFAnalysis_FDR0.10/glmQLF_2WayANOVA_TvsN_topTags_filtered.csv")
SETInterIn <- geneCountsInter[,1]
SETTreatIn <- geneCountsTreat[,1]
SETTolIn <- geneCountsTol[,1]

SETInter <- SETInterIn[SETInterIn %in% SETDDR]
SETTreat <- SETTolIn[SETTolIn %in% SETDDR]
SETTol <- SETTolIn[SETTolIn %in% SETDDR]

lt = list(Interaction=SETInter, Treatment=SETTreat, Tolerance=SETTol)
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

lt = list(Interaction=SETInter, Treatment=SETTreat, Tolerance=SETTol, Tribolium=SETTCAST, Drosophila=SETDMEL)
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
