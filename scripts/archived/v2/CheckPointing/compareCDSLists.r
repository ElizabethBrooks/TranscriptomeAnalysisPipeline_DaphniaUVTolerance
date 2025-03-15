setwd("/Users/bamflappy/PfrenderLab/PA42_v4.1")

olympicsList <- read.csv(file="Olympics_longest_cds_list.txt", header=FALSE)
paList <- read.csv(file="PA42_v4.1_longest_cds_list.txt", header=FALSE)

#Compare the two sets of longest CDS
length(setdiff(olympicsList[,1], paList[,1]))
length(setdiff(paList[,1], olympicsList[,1]))
length(intersect(paList[,1], olympicsList[,1]))
