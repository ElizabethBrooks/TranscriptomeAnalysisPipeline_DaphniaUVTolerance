# script to compare WGCN for a condition

# install libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MODA")

# load libraries
library("MODA")

# import expression data
countsTable <- read.csv(file="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/GeneCountsAnalyzed_KAP4/geneCounts_merged_genome_counted_htseq_run1.txt", row.names="gene", header=TRUE)[ ,1:24]

# subset expression data for expression profile 1
exp_profile1

# subset expression data for expression profile 2
exp_profile2

# initialize variables for input to MODA functions
WGCNResultFolder = 'OLYM_WGCN' # where middle files are stored
CuttingCriterion = 'Density' # could be Density or Modularity
indicator_profile1 = 'All'     # indicator for data profile 1
indicator_profile2 = 'UV'      # indicator for data profile 2
specificTheta = 0.1 #threshold to define condition specific modules
conservedTheta = 0.1#threshold to define conserved modules

# detect modules for network 1
intModules_profile1 <- WeightedModulePartitionHierarchical(exp_profile1,
                                                           WGCNResultFolder,
                                                           indicator_profile1,
                                                           CuttingCriterion)

# detect modules for network 2
intModules_profile2 <- WeightedModulePartitionHierarchical(exp_profile2,
                                                           WGCNResultFolder,
                                                           indicator_profile2,
                                                           CuttingCriterion)

# network comparison
CompareAllNets(WGCNResultFolder,
               intModules_profile1,
               indicator_profile1,
               intModules_profile2,
               indicator_profile2,
               specificTheta,
               conservedTheta)