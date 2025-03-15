#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install()

#Load the edgeR library
library(edgeR)
library(GO.db)
library(reshape2)
library(ggplot2)

#Import gene count data
countsTable <- read.csv(file= "daphniaRawCounts_onlyEntrez_noDups.csv", 
                        row.names="gene")[ ,31:36]
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))
#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table

#Perform GO analysis
go <- goana(tested, geneid = rownames(tested))
top_go <- topGO(go) #add sort = 'up' ? if you want only the GO terms that are overrepresented

#list of significant over/under represented GO terms for each genotype
go_GENOTYPE.up <- go[go$P.Up <= 0.05, ]
go_GENOTYPE.down <- go[go$P.Down <= 0.05, ]

#See intersection of GO terms
tolerant_intersect.up <- intersect(go_Y023.up$Term, go_R2.up$Term)
tolerant_intersect.down <- intersect(go_Y023.down$Term, go_R2.down$Term)
not_tolerant_intersect.up <- intersect(go_E05.up$Term, go_Y05.up$Term)
not_tolerant_intersect.down <- intersect(go_E05.down$Term, go_Y05.down$Term)
intersect_all.up <- intersect(intersect(go_Y023.up$Term, go_R2.up$Term), 
                              intersect(go_E05.up$Term, go_Y05.up$Term))
intersect_all.down <- intersect(intersect(go_E05.down$Term, go_Y05.down$Term), 
                                intersect(go_Y023.down$Term, go_R2.down$Term))

#Perform KEGG analysis
keg <- kegga(tested, geneid = rownames(tested))
top_keg <- topKEGG(keg)


#Plot % of genes in top GO term that were DE'd 
top_go_plot <- data.frame(top_go$Term, (100* top_go$Up/top_go$N), (100 * top_go$Down/top_go$N))
colnames(top_go_plot) <- c('GOTerm', 'Up', 'Down')

d = melt(top_go_plot, id.vars = "GOTerm")

jpeg(file = 'GO_plot_GENOTYPE.jpg', width = 800, height = 600)
p <- ggplot(data = d,
       mapping = aes(x = GOTerm, y = value, fill = variable)) + 
  geom_col(position = position_dodge()) 
p <- p + theme(axis.text.x = element_text(angle = 90)) 
p <- p + ylab('Percent of genes') + xlab('GO Term') + labs(fill = 'Regulation')
p
dev.off()



