#install.packages('VennDiagram')
library(VennDiagram)

#Import data
ids <- read.csv("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/mergedTopGenesTags_exactTest.csv")
#Define sets for diagram
SET1 <- ids$Y05
SET2 <- ids$Y023
SET3 <- ids$E05
SET4 <- ids$R2
SET5 <- ids$PA
SET6 <- ids$Sierra

#Draw the diagram from the Olympics sets
v1 <- venn.diagram(list(Y05=SET1, Y023=SET2, E05=SET3, R2=SET4),
                  fill = c("red", "green", "white", "blue"),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/plotOlympicsVennn_topTags.jpg")
grid.newpage()
grid.draw(v1)
dev.off()

#Draw the diagram comparing the Olympics and PA sets
v2 <- venn.diagram(list(Y05=SET1, Y023=SET2, E05=SET3, R2=SET4, PA=SET5),
                  fill = c("red", "green", "white", "blue", "yellow"),
                  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/plotOlympicsPAVennn_topTags.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

#Draw the diagram comparing the Olympics and Sierra sets
v3 <- venn.diagram(list(Y05=SET1, Y023=SET2, E05=SET3, R2=SET4, Sierra=SET6),
                   fill = c("red", "green", "white", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/plotOlympicsSierraVennn_topTags.jpg")
grid.newpage()
grid.draw(v3)
dev.off()

#Draw the diagram comparing the tolerant sets
v4 <- venn.diagram(list(Y05=SET1, E05=SET3, Sierra=SET6),
                   fill = c("red", "green","blue"),
                   alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/plotTolerantVenn_topTags.jpg")
grid.newpage()
grid.draw(v4)
dev.off()

#Draw the diagram comparing the non-tolerant sets
v5 <- venn.diagram(list(Y023=SET2, R2=SET4, PA=SET5),
                   fill = c("red", "green","blue"),
                   alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GeneCounts_Stats/plotNonTolerantVennn_topTags.jpg")
grid.newpage()
grid.draw(v5)
dev.off()