# script to compare WGCN for a condition

#Set working directory
workingDir="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4"
setwd(workingDir); 

# install libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MODA")

# load libraries
library("MODA")

# where middle files are stored
WGCNResultFolder <- "OLYM_WGCN_MODA"

# could be Density or Modularity
CuttingCriterion <- "Density"

# threshold to define condition specific modules
specificTheta <- 0.1

# threshold to define conserved modules
conservedTheta <- 0.1

# setup indicators for each expression profile
ind_all <- "All"
# treatment
ind_UV <- "UV"
ind_VIS <- "VIS"
# tolerance
ind_tol <- "tol"
ind_nTol <- "nTol"
# genotype
ind_Y05 <- "Y05"
ind_Y023 <- "Y023"
ind_R2 <- "R2"
ind_E05 <- "E05"

# import expression data
countsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/DEGenotypes/glmQLF_normalizedCounts.csv", row.names="gene", header=TRUE)[ ,1:24]
# treamtent
exp_UV <- countsTable[,grepl("UV", names(countsTable))]
exp_VIS <- countsTable[,grepl("VIS", names(countsTable))]
# tolerance
exp_tol <- countsTable[,grepl("Y05|Y023", names(countsTable))]
exp_nTol <- countsTable[,grepl("R2|E05", names(countsTable))]
# genotype
exp_Y05 <- countsTable[,grepl("Y05", names(countsTable))]
exp_Y023 <- countsTable[,grepl("Y023", names(countsTable))]
exp_R2 <- countsTable[,grepl("R2", names(countsTable))]
exp_E05 <- countsTable[,grepl("E05", names(countsTable))]

# view correlation matrices for each network
heatmap(cor(as.matrix(countsTable)))
# treatment
heatmap(cor(as.matrix(exp_UV)))
heatmap(cor(as.matrix(exp_VIS)))
# tolerance
heatmap(cor(as.matrix(exp_tol)))
heatmap(cor(as.matrix(exp_nTol)))
# genotype
heatmap(cor(as.matrix(exp_Y05)))
heatmap(cor(as.matrix(exp_Y023)))
heatmap(cor(as.matrix(exp_R2)))
heatmap(cor(as.matrix(exp_E05)))

# detect modules for the network of each expression profile
intModules_all <- WeightedModulePartitionHierarchical(countsTable,
                                                           WGCNResultFolder,
                                                           ind_all,
                                                           CuttingCriterion)
# treatment
intModules_UV <- WeightedModulePartitionHierarchical(exp_UV,
                                                           WGCNResultFolder,
                                                           ind_UV,
                                                           CuttingCriterion)
intModules_VIS <- WeightedModulePartitionHierarchical(exp_VIS,
                                                     WGCNResultFolder,
                                                     ind_VIS,
                                                     CuttingCriterion)
# tolerance
intModules_tol <- WeightedModulePartitionHierarchical(exp_tol,
                                                     WGCNResultFolder,
                                                     ind_tol,
                                                     CuttingCriterion)
intModules_nTol <- WeightedModulePartitionHierarchical(exp_nTol,
                                                     WGCNResultFolder,
                                                     ind_nTol,
                                                     CuttingCriterion)
# genotype
intModules_Y05 <- WeightedModulePartitionHierarchical(exp_Y05,
                                                     WGCNResultFolder,
                                                     ind_Y05,
                                                     CuttingCriterion)
intModules_Y023 <- WeightedModulePartitionHierarchical(exp_Y023,
                                                     WGCNResultFolder,
                                                     ind_Y023,
                                                     CuttingCriterion)
intModules_R2 <- WeightedModulePartitionHierarchical(exp_R2,
                                                     WGCNResultFolder,
                                                     ind_R2,
                                                     CuttingCriterion)
intModules_E05 <- WeightedModulePartitionHierarchical(exp_E05,
                                                     WGCNResultFolder,
                                                     ind_E05,
                                                     CuttingCriterion)

# network comparison for each subset to the entire background set
# treatment
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_UV,
               ind_UV,
               specificTheta,
               conservedTheta)
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_VIS,
               ind_VIS,
               specificTheta,
               conservedTheta)
# tolerance
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_tol,
               ind_tol,
               specificTheta,
               conservedTheta)
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_nTol,
               ind_nTol,
               specificTheta,
               conservedTheta)
# genotype
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_Y05,
               ind_Y05,
               specificTheta,
               conservedTheta)
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_Y023,
               ind_Y023,
               specificTheta,
               conservedTheta)
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_R2,
               ind_R2,
               specificTheta,
               conservedTheta)
CompareAllNets(WGCNResultFolder,
               intModules_all,
               ind_all,
               intModules_E05,
               ind_E05,
               specificTheta,
               conservedTheta)
