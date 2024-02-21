#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install("gridExtra")

#Load libraries
library(ggplot2)
library(gridExtra)
library(rcartocolor)
library(dplyr)

# Plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[5], plotColors[6])

# turn off scientific notation
options(scipen = 999)


# DE data
# retrieve working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis"

# set working directory
setwd(workingDir)

# read in data on significant GO terms (BP, MF, and CC) for each effect
treatment_BP_GO_terms <- read.csv('treatment_BP_sigGO_terms.csv')
#treatment_MF_GO_terms <- read.csv('treatment_MF_sigGO_terms.csv')
#treatment_CC_GO_terms <- read.csv('treatment_CC_sigGO_terms.csv')

## treatment effect
# filter for top 5 significant terms
treatment_BP_GO_sig <- treatment_BP_GO_terms#[1:5, ]
#treatment_MF_GO_sig <- treatment_MF_GO_terms#[1:5, ]
#treatment_CC_GO_sig <- treatment_CC_GO_terms#[1:5, ]

# add a column labeling the effect to each GO term set
treatment_BP_plot_table <- cbind("Set" = 'Treatment', treatment_BP_GO_sig)
#treatment_MF_plot_table <- cbind("Set" = 'Treatment', treatment_MF_GO_sig)
#treatment_CC_plot_table <- cbind("Set" = 'Treatment', treatment_CC_GO_sig)

# add GO level tags
treatment_BP_plot_table <- cbind('Level' = 'BP', treatment_BP_plot_table)
#treatment_MF_plot_table <- cbind('Level' = 'MF', treatment_MF_plot_table)
#treatment_CC_plot_table <- cbind('Level' = 'CC', treatment_CC_plot_table)

# plotting
#all_plot_table <- rbind(treatment_BP_plot_table, treatment_MF_plot_table, treatment_CC_plot_table)

# select columns
plot_tableDE <- select(treatment_BP_plot_table, Level, Set, GO.ID, Term, Significant, weightFisher)


# module data
# retrieve working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30"

# set working directory
setwd(workingDir)

# retrieve subset tag
set <- "OLYM"

# set the minimum module size
minModSize <- "30"

# retrieve WGCNA directory
inDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"

# set the full subset tag name
tag <- paste(set, minModSize, sep="_")

# Load network data saved in the second part
importFile <- paste(tag, "networkConstruction-stepByStep.RData", sep="-")
importFile <- paste(inDir, importFile, sep="/")
lnames2 = load(file = importFile)

# create list of module colors mapped to numbers
num_mods <- length(unique(moduleColors)) + 1
color_list <- c(unique(moduleColors), "None")

# create data frame to hold plotting data for each module
module_BP_results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(module_BP_results) <- c("Set","Term","Significant","weightFisher")
#module_MF_results = data.frame(matrix(ncol = 4, nrow = 0))
#colnames(module_MF_results) <- c("Set","Term","Significant","weightFisher")
#module_CC_results = data.frame(matrix(ncol = 4, nrow = 0))
#colnames(module_CC_results) <- c("Set","Term","Significant","weightFisher")

# loop through each module color
for(j in 1:num_mods){
  # read in data on significant GO terms (BP, MF, and CC) 
  BP_GO_sig <- read.csv(paste(color_list[j], 'BP_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  #MF_GO_sig <- read.csv(paste(color_list[j], 'MF_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  #CC_GO_sig <- read.csv(paste(color_list[j], 'CC_sigGO_terms.csv', sep="_"), row.names = NULL)#[1:5,]
  
  # add a column labeling the color to each GO term set
  BP_GO_sig <- cbind("Set" = color_list[j], BP_GO_sig[,1], BP_GO_sig[,2], BP_GO_sig[,4], BP_GO_sig[,6])
  #MF_GO_sig <- cbind("Set" = color_list[j], MF_GO_sig[,1], MF_GO_sig[,2], MF_GO_sig[,4], MF_GO_sig[,6])
  #CC_GO_sig <- cbind("Set" = color_list[j], CC_GO_sig[,1], CC_GO_sig[,2], CC_GO_sig[,4], CC_GO_sig[,6])
  
  # combine all tables
  module_BP_results <- rbind(module_BP_results, BP_GO_sig)
  #module_MF_results <- rbind(module_MF_results, MF_GO_sig)
  #module_CC_results <- rbind(module_CC_results, CC_GO_sig)
}

# re-name columns
colnames(module_BP_results) <- c("Set","GO.ID","Term","Significant","weightFisher")
#colnames(module_MF_results) <- c("Set","GO.ID","Term","Significant","weightFisher")
#colnames(module_CC_results) <- c("Set","GO.ID","Term","Significant","weightFisher")

# add GO level tags
module_BP_results <- cbind('Level' = 'BP', module_BP_results)
#module_MF_results <- cbind('Level' = 'MF', module_MF_results)
#module_CC_results <- cbind('Level' = 'CC', module_CC_results)

# combine GO level results for plotting
#plotTableModule <- rbind(module_BP_results, module_MF_results, module_CC_results)

# remove NAs
plotTableModule <- na.omit(module_BP_results)

# remove < signs
plotTableModule$weightFisher <- gsub("<","",as.character(plotTableModule$weightFisher))

# select columns
plot_tableModule <- select(plotTableModule, Level, Set, GO.ID, Term, Significant, weightFisher)

# combine tables
plot_table <- rbind(plot_tableDE, plot_tableModule)

# geneIDs with repair or other interesting GO terms
# The repair pathway GO terms we specifically explored in the analysis were DNA repair (GO:0006281), mismatch repair (MMR; GO:0006298), base excision repair (BER; GO:0006284), homologous recombination (HR; GO:0035825), nucleotide excision repair (NER; GO:0006289), intrastrand crosslink repair (ICL repair; GO:0036297), double strand break repair (DSBR; GO:0006302), single strand break repair (SSBR; GO:0000012).
repairTerms <- list(DNAR = "GO:0006281",
                    MMR = "GO:0006298",
                    BER = "GO:0006284",
                    HR = "GO:0035825",
                    NER = "GO:0006289",
                    ICLR = "GO:0036297",
                    DSBR = "GO:0006302",
                    SSBR = "GO:0000012"
)

# The radiation response terms included response to radiation (GO:0009314), cellular response to radiation (GO:0071478), phototransduction UV (GO:0007604),
# response to UV (GO:0009411), response to UV-A (GO:0070141), detection of UV (GO:0009589), cellular response to UV (GO:0034644), cellular response to UV-A (GO:0071492), 
# response to UV-B4 (GO:0010224), cellular response to UV-B (GO:0071493), response to UV-C (GO:0010225), regulation of mRNA stability involved in cellular response to UV (GO:1902629), 
# cellular response to UV-C (GO:0071494), regulation of translation involved in cellular response to UV (GO:1904803).
radiationTerms <- list(radiation = "GO:0009314",
                       ionizing = "GO:0010212",
                       gamma = "GO:0010332",
                       cell = "GO:0071478",
                       cellGamma = "GO:0071480",
                       cellIonizing = "GO:0071479",
                       regGamma = "GO:2001228",
                       regCellGamma = "GO:1905843",
                       xRay = "GO:0010165",
                       cellXRay = "GO:0071481",
                       regCellXRay = "GO:2000683",
                       photoUV = "GO:0007604",
                       UV = "GO:0009411",
                       UVA = "GO:0070141",
                       detectUV = "GO:0009589",
                       cellUV = "GO:0034644",
                       cellUVA = "GO:0071492",
                       UVB4 = "GO:0010224",
                       cellUVB = "GO:0071493",
                       UVC = "GO:0010225",
                       mRNACellUV = "GO:1902629",
                       cellUVC = "GO:0071494",
                       transCellUV = "GO:1904803"
)

# The stress terms included response to stress (GO:0006950), response to oxidative stress (GO:0006979), cellular response to oxidative stress (GO:0034599), 
# cellular response to reactive oxygen (GO:0034614), regulation of translation in response to oxidative stress (GO:0043556), regulation of cellular response to oxidative stress (GO:1900407).
stressTerms <- list(stress = "GO:0006950",
                    oxidative = "GO:0006979",
                    cellOxidative = "GO:0034599",
                    cellReactiveOxy = "GO:0034614",
                    regTransOxidative = "GO:0043556",
                    MAPKKK = "GO:1990315",
                    hydroperoxide = "GO:0071447",
                    senescence = "GO:0090403",
                    symbiont = "GO:0052164",
                    regCellOxidative = "GO:1900407"
)

# combine candidate GO term lists
candidateTerms <- c(repairTerms, radiationTerms, stressTerms)

# initilize data frame
plot_tableSubset <- plot_table[plot_table$GO.ID %in% candidateTerms,]

# select the salmon4 and skyblue modules and Treatment effect
plot_tableSubset <- plot_tableSubset[plot_tableSubset$Set %in% c("Treatment","salmon4","skyblue"),]

# setup DE facet groups
x_axis_order <- factor(plot_tableSubset$Set, levels = c('Treatment', 'salmon4', 'skyblue'))
#facet <- factor(plot_tableSubset$Level, levels = c('BP', 'CC', 'MF'))

# create dot plot of significant GO terms
dotplot <- ggplot(data = plot_tableSubset, aes(x = x_axis_order, y = GO.ID, size = as.numeric(Significant), color = as.numeric(weightFisher))) + 
  #facet_grid(rows = facet, space = 'free_y', scales = 'free') +
  geom_point() +
  #scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
  scale_color_gradientn(colors = plotColorSubset) +
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, color = plot_tableSubset$selection)) +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_text(position = position_dodge(width = 1), aes(x=effect, y=0)) +
  xlab('Color') +
  ylab('GO Term') + 
  labs(color = 'P-Value', size = 'Gene Rank')

# view plot
dotplot

# save the plot to a PDF file
#ggsave('dotplotDEModule_sigGO.pdf', plot = dotplot, device = 'pdf')
