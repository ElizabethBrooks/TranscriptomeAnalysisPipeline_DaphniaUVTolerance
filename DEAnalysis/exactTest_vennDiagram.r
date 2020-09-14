##TODO: Add percents

#install.packages('VennDiagram')
library(VennDiagram)

#Import data
counts <- read.csv(file=args[1])
#Define sets for diagram
SET1 <- counts$Y05
SET2 <- counts$Y023
SET3 <- counts$E05
SET4 <- counts$R2
SET5 <- counts$PA
SET6 <- counts$Sierra

#Draw the diagram from the Olympics sets
v1 <- venn.diagram(list(Y05=SET1, E05=SET3, Y023=SET2, R2=SET4),
                  fill = c("red", "green", "white", "blue"),
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
jpeg("plotOlympicsVennn_exactTest.jpg")
grid.newpage()
grid.draw(v1)
dev.off()

#Draw the diagram comparing the Olympics and PA sets
v2 <- venn.diagram(list(Y05=SET1, E05=SET3, Y023=SET2, R2=SET4, PA=SET5),
                  fill = c("red", "green", "white", "blue", "yellow"),
                  alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),
                  filename=NULL)
jpeg("plotOlympicsPAVennn_exactTest.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

#Draw the diagram comparing the Olympics and Sierra sets
v3 <- venn.diagram(list(Y05=SET1, E05=SET3, Y023=SET2, R2=SET4, Sierra=SET6),
                   fill = c("red", "green", "white", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("plotOlympicsSierraVennn_exactTest.jpg")
grid.newpage()
grid.draw(v3)
dev.off()

#Draw the diagram comparing the tolerant sets
v4 <- venn.diagram(list(Y05=SET1, E05=SET3, Sierra=SET6),
                   fill = c("red", "green","blue"),
                   alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("plotTolerantVenn_exactTest.jpg")
grid.newpage()
grid.draw(v4)
dev.off()

#Draw the diagram comparing the non-tolerant sets
v5 <- venn.diagram(list(Y023=SET2, R2=SET4, PA=SET5),
                   fill = c("red", "green","blue"),
                   alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("plotNonTolerantVennn_exactTest.jpg")
grid.newpage()
grid.draw(v5)
dev.off()

#Draw the diagram comparing the Olympic tolerant sets
v6 <- venn.diagram(list(Y05=SET1, E05=SET3),
                   fill = c("red","blue"),
                   alpha = c(0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("plotOlympicTolerantVenn_exactTest.jpg")
grid.newpage()
grid.draw(v6)
dev.off()

#Draw the diagram comparing the Olympic non-tolerant sets
v7 <- venn.diagram(list(Y023=SET2, R2=SET4),
                   fill = c("red","blue"),
                   alpha = c(0.5, 0.5), cat.cex = 1.5, cex=1.5,
                   filename=NULL)
jpeg("plotOlympicNonTolerantVennn_exactTest.jpg")
grid.newpage()
grid.draw(v7)
dev.off()
