#install.packages('VennDiagram')
library(VennDiagram)

#Import data
tags <- read.csv(file="/home/mae/Documents/RNASeq_Workshop_ND/GLMQLFvsGSEA/vennSet_GLMvsGSEA.csv")
#Define sets for diagram
SET1 <- tags$GLM_Treat
SET2 <- tags$GLM_Tol
SET3 <- tags$GLM_Inter
SET4 <- tags$GSEA
#Replace NAs
SET1[is.na(SET1)] <- ""
SET2[is.na(SET2)] <- ""
SET3[is.na(SET3)] <- ""
SET4[is.na(SET4)] <- ""

#Draw the diagram from the Olympics sets
v1 <- venn.diagram(list(GLM.Treat=SET1, GLM.Tol=SET2, GLM.Inter=SET3, GSEA=SET4),
                   fill = c("red", "green", "blue", "white"),
                   alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GLMQLFvsGSEA/plotResults_GLMvsGSEA.jpg")
grid.newpage()
grid.draw(v1)
dev.off()

#Draw the diagram comparing the tolerant sets
v2 <- venn.diagram(list(GLM.Treat=SET1, GLM.Tol=SET2, GLM.Inter=SET3),
                   fill = c("red", "green", "blue"),
                   alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("/home/mae/Documents/RNASeq_Workshop_ND/GLMQLFvsGSEA/plotResults_GLM.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

