#!/usr/bin/env Rscript
#Usage: Rscript geneSetTest_camera.r countsFile factorGroupingFile
#Usage Ex: Rscript geneSetTest_camera.r PA42_v4.1_normalizedCountsOlympics_uniprot.csv expDesign_camera_Olympics.csv
#R script to perform gene set enrichment testing using camera

# Specify a directory to save the output
output.directory <- "/Users/bryanmichalek/Documents/Notre_Dame/Spring 2022/Pfrender/CAMERA"

# Load libraries
library("limma")
library("msigdbr")
library("dplyr")
library("tidyverse")
library("ggnewscale")

#Import gene count data
#Row names are set to the uniprot IDs
# to match the gene IDs of the MSigDB KEGG DNA repair gene sets
#countsTable <- read.csv(file=args[1])
countsTable <- read.csv(file="PA42_v4.1_normalizedLogCountsOlympics_uniprot.csv")

#Create a subset of the input counts table containing only the gene counts
counts <- countsTable[3:26]

#Import grouping factor
#targets <- read.csv(file=args[2], row.names="sample")
targets <- read.csv(file="expDesign_camera_Olympics.csv", row.names="sample")

#Setup a design matrix
tolerance <- targets$tolerance
treatment <- targets$treatment

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + tolerance + treatment + tolerance:treatment)
colnames(design) <- c("(Intercept)","NTol.Tol","UV.VIS","Interaction")

#--------------------------------------------------------------------------------------------------------------------------------
#Import gene sets
#Here, I use the msigdbr function to grab all the KEGG pathways
#Then, store the names of gene sets in a vector which will be used as the names in the list() function
#The values in the list will be the indices of countsTable that are contained in each respective gene set name. 
#This is accomplished with the for loop (which gets the indices and then assigns them to the correct name within the list...iterating over ever possible name)
CP.KEGG <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"))
GO.BP <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
GO.CC <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"))
GO.MF <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"))

geneSetNames.CP.KEGG <- unique(CP.KEGG$gs_name)
geneSetNames.GO.BP <- unique(GO.BP$gs_name)
geneSetNames.GO.CC <- unique(GO.CC$gs_name)
geneSetNames.GO.MF <- unique(GO.MF$gs_name)

set_geneSet_indices <- function(geneList, geneSetsNames){
  geneSetList <- list()
  for(set.name in geneSetsNames){
    index <- which(countsTable$sprot %in% geneList[geneList$gs_name == set.name, "gene_symbol"])
    geneSetList[[set.name]] <- index
  }
  return(geneSetList)
}

geneSets.CP.KEGG <- set_geneSet_indices(CP.KEGG, geneSetNames.CP.KEGG)
geneSets.GO.BP <- set_geneSet_indices(GO.BP, geneSetNames.GO.BP)
geneSets.GO.CC <- set_geneSet_indices(GO.CC, geneSetNames.GO.CC)
geneSets.GO.MF <- set_geneSet_indices(GO.MF, geneSetNames.GO.MF)

#--------------------------------------------------------------------------------------------------------------------------------
#Summary table containing counts of genes in each set (will have 2 columns):
#The number of total genes in each set from msigdbr (regardless of whether each is found in our countsTable)
count_genes_long <- function(geneList){
  count_long <- count(geneList, gs_name)
  colnames(count_long) <- c("gs_name", "Total")
  return(count_long)
}

CP.KEGG.total.count <- count_genes_long(CP.KEGG)
GO.BP.total.count <- count_genes_long(GO.BP)
GO.CC.total.count <- count_genes_long(GO.CC)
GO.MF.total.count <- count_genes_long(GO.MF)

#The number of genes in our countsTable that were found in the complete list for each gene set
count_genes_short <- function(geneList){
  count_short <- lapply(geneList, length) %>% unlist() 
  count_short <- data.frame(gs_name = names(count_short), n = count_short, row.names = 1:length(count_short))
  colnames(count_short) <- c("gs_name", "Matches")
  return(count_short)
}

CP.KEGG.matched.count <- count_genes_short(geneSets.CP.KEGG)
GO.BP.matched.count <- count_genes_short(geneSets.GO.BP)
GO.CC.matched.count <- count_genes_short(geneSets.GO.CC)
GO.MF.matched.count <- count_genes_short(geneSets.GO.MF)

#Merge total counts with the count of matches for each database
CP.KEGG.counts.final <- merge(CP.KEGG.total.count, CP.KEGG.matched.count, by = "gs_name") %>% column_to_rownames(var = "gs_name")
GO.BP.counts.final <- merge(GO.BP.total.count, GO.BP.matched.count, by = "gs_name") %>% column_to_rownames(var = "gs_name")
GO.CC.counts.final <- merge(GO.CC.total.count, GO.CC.matched.count, by = "gs_name") %>% column_to_rownames(var = "gs_name")
GO.MF.counts.final <- merge(GO.MF.total.count, GO.MF.matched.count, by = "gs_name") %>% column_to_rownames(var = "gs_name")

#--------------------------------------------------------------------------------------------------------------------------------
#Batch gene set tests
#Inter-gene correlation will be estimated for each tested set
cam <- function(geneSets){
  #Tolerance contrast
  df_1 <- camera(counts, geneSets, design, contrast = 2, inter.gene.cor = NA)
  #Treatment contrast
  df_2 <- camera(counts, geneSets, design, contrast = 3, inter.gene.cor = NA)
  df_2 <- df_2[, -1:-2]
  #Interaction contrast
  df_3 <- camera(counts, geneSets, design, contrast = 4, inter.gene.cor = NA)
  df_3 <- df_3[, -1:-2]
  #Combine and return results
  df_final <- merge(df_1, df_2, by = "row.names", suffixes = c("_toler", "_treat")) %>%
    merge(df_3, by.x = "Row.names", by.y = 0)
  colnames(df_final)[c(1, 10:12)] <- c("Gene_set", "Direction_inter", "PValue_inter", "FDR_inter")
  return(df_final)
}
CP.KEGG.cam <- cam(geneSets.CP.KEGG)
GO.BP.cam <- cam(geneSets.GO.BP)
GO.CC.cam <- cam(geneSets.GO.CC)
GO.MF.cam <- cam(geneSets.GO.MF)

#-------------------------------------------------------------------------------------
#Pre-ranked (PR) gene set test version
#Fit the linear model using the design matrix
fit <- eBayes(lmFit(counts, design))
#Function for running pre-ranked tests
camPR <- function(geneSets){
  #Using moderated t-statistic for tolerance effect
  df_1 <- cameraPR(fit$t[,2], geneSets, use.ranks = TRUE)
  #Using moderated t-statistic for treatment effect
  df_2 <- cameraPR(fit$t[,3], geneSets, use.ranks = TRUE)[, -1] #Removing NGenes column since it is the same as in df_1
  #Using moderated t-statistic for interaction effect
  df_3 <- cameraPR(fit$t[,4], geneSets, use.ranks = TRUE)[, -1] #Removing NGenes column since it is the same as in df_1
  #Using moderated F-statistic
  df_4 <- cameraPR(fit$F, geneSets, use.ranks = TRUE)[, -1] #Removing NGenes column since it is the same as in df_1
  #Combine and return
  df_final <- merge(df_1, df_2, by = "row.names", suffixes = c("_toler", "_treat")) %>% 
    merge(df_3, by.x = "Row.names", by.y = 0) %>% 
    merge(df_4, by.x = "Row.names", by.y = 0, suffixes = c("_inter", "_F"))
  colnames(df_final)[1] <- "Gene_set"
  return(df_final)
}

CP.KEGG.camPR <- camPR(geneSets.CP.KEGG)
GO.BP.camPR <- camPR(geneSets.GO.BP)
GO.CC.camPR <- camPR(geneSets.GO.CC)
GO.MF.camPR <- camPR(geneSets.GO.MF)

#-------------------------------------------------------------------------------------
# geneSetTest for non-directional tests (use sapply with FUN = geneSetTest and then list all the optional arguments to geneSetTest after this)
gst <- function(geneSets){
  # Using moderated t-statistic for tolerance
  df_1 <- sapply(geneSets, FUN = geneSetTest, statistics = fit$t[, 2], alternative = "mixed", type = "t") %>% as.data.frame()
  colnames(df_1) <- "PValue_toler"
  # Using moderated t-statistic for treatment
  df_2 <- sapply(geneSets, FUN = geneSetTest, statistics = fit$t[, 3], alternative = "mixed", type = "t") %>% as.data.frame()
  colnames(df_2) <- "PValue_treat"
  # Using moderated t-statistic for interaction
  df_3 <- sapply(geneSets, FUN = geneSetTest, statistics = fit$t[, 4], alternative = "mixed", type = "t") %>% as.data.frame()
  colnames(df_3) <- "PValue_inter"
  # Using moderated F-statistic
  df_4 <- sapply(geneSets, FUN = geneSetTest, statistics = fit$F, alternative = "mixed", type = "f") %>% as.data.frame()
  colnames(df_4) <- "PValue_F"
  # Combine and return
  df_final <- merge(df_1, df_2, by = "row.names") %>% 
    merge(df_3, by.x = "Row.names", by.y = 0) %>% 
    merge(df_4, by.x = "Row.names", by.y = 0)
  colnames(df_final)[1] <- "Gene_set"
  return(df_final)
}
CP.KEGG.gst <- gst(geneSets.CP.KEGG)
GO.BP.gst <- gst(geneSets.GO.BP)
GO.CC.gst <- gst(geneSets.GO.CC)
GO.MF.gst <- gst(geneSets.GO.MF)

#--------------------------------------------------------------------------------------------------------------------------------
# Process the 12 data frames containing results of each test to prepare them for plotting
# Write function to format data for the plot
reformat_results <- function(database, analysis.function, df_results){
  # Add column for what database the geneSets belong to
  df_results["Database"] <- database
  df_results <- df_results %>% relocate(Database)
  
  # cam()/camPR() have different data frame outputs from gst(). Need to process them slightly different. 
  if(analysis.function == "gst"){
    df_final <- df_results %>% pivot_longer(
      cols = c(PValue_toler, PValue_treat, PValue_inter, PValue_F),
      names_to = "ANOVA_group",
      names_prefix = "PValue.",
      values_to = "PValue"
    )
  }else if(analysis.function %in% c("cam", "camPR")){
    if(analysis.function == "cam"){
      columns <- c("Direction_toler", "PValue_toler", "FDR_toler",
                   "Direction_treat", "PValue_treat", "FDR_treat",
                   "Direction_inter", "PValue_inter", "FDR_inter")
    }else{
      columns <- c("Direction_toler", "PValue_toler", "FDR_toler",
                   "Direction_treat", "PValue_treat", "FDR_treat",
                   "Direction_inter", "PValue_inter", "FDR_inter",
                   "Direction_F", "PValue_F", "FDR_F")
    }
    
    df_final <- df_results %>% pivot_longer(
      cols = columns,
      names_to = all_of(c("statistic", "ANOVA_group")),
      names_sep = "_",
      values_to = "values",
      values_transform = list(values = as.character)
    ) %>% pivot_wider(
      names_from = "statistic",
      values_from = "values"
    )
    df_final$PValue <- as.numeric(df_final$PValue)
    df_final$FDR <- as.numeric(df_final$FDR)
    df_final <- df_final %>% select(-NGenes)
  }
  return(df_final)
}

  # Prepare inputs to run the function many times
analysis.functions <- rep(c("cam", "camPR", "gst"), each = 4)
databases <- rep(c("KEGG", "BP", "CC", "MF"), 3)
  # Create list of all of our data frames with results
results.list <- list(CP.KEGG.cam, GO.BP.cam, GO.CC.cam, GO.MF.cam,
                     CP.KEGG.camPR, GO.BP.camPR, GO.CC.camPR, GO.MF.camPR,
                     CP.KEGG.gst, GO.BP.gst, GO.CC.gst, GO.MF.gst)
  # List of the transformed data frames we will use to plot
results.list.reformatted <- vector(mode = "list", length = 12)
names(results.list.reformatted) <- c("CP.KEGG.cam.reformatted", "GO.BP.cam.reformatted", "GO.CC.cam.reformatted", "GO.MF.cam.reformatted",
                      "CP.KEGG.camPR.reformatted", "GO.BP.camPR.reformatted", "GO.CC.camPR.reformatted", "GO.MF.camPR.reformatted",
                      "CP.KEGG.gst.reformatted", "GO.BP.gst.reformatted", "GO.CC.gst.reformatted", "GO.MF.gst.reformatted")
  # Run the function 12 times to transform each data frame
for(i in 1:12){
  results.list.reformatted[[i]] <- reformat_results(databases[i], analysis.functions[i], results.list[[i]])
  
  # Add the counts information from above to each data frame so we can use these columns as the "size" factor on plots
  if(databases[i] == "KEGG"){
    results.list.reformatted[[i]] <- merge(results.list.reformatted[[i]], CP.KEGG.counts.final, by.x = "Gene_set", by.y = "row.names", all.x = TRUE)
  }else if(databases[i] == "BP"){
    results.list.reformatted[[i]] <- merge(results.list.reformatted[[i]], GO.BP.counts.final, by.x = "Gene_set", by.y = "row.names", all.x = TRUE)
  }else if(databases[i] == "CC"){
    results.list.reformatted[[i]] <- merge(results.list.reformatted[[i]], GO.CC.counts.final, by.x = "Gene_set", by.y = "row.names", all.x = TRUE)
  }else{
    results.list.reformatted[[i]] <- merge(results.list.reformatted[[i]], GO.MF.counts.final, by.x = "Gene_set", by.y = "row.names", all.x = TRUE)
  }
}

#--------------------------------------------------------------------------------------------------------------------------------
# Create the tables we need to plot. There will be 3 of them (cam, camPR, and gst)
# We will create our list of gene sets to display on the y-axis by grabbing the 4 most significant gene sets in each of 
# the treatment, tolerant, interaction, and F stat groups for each of the databases (KEGG, BP, CC, MF)
create_plot_table <- function(results.full, analysis.function, n.significant){
  results.subset <- results.full[which(analysis.functions == analysis.function)]
  unique.sets <- c()
  plot.table <- data.frame()
  
  for(i in 1:length(results.subset)){
    treat <- results.subset[[i]] %>% filter(ANOVA_group == "treat") %>% 
      mutate(Rank = min_rank(PValue)) %>% 
      filter(Rank <= n.significant) %>% 
      select(Gene_set)
    
    toler <- results.subset[[i]] %>% filter(ANOVA_group == "toler") %>% 
      mutate(Rank = min_rank(PValue)) %>% 
      filter(Rank <= n.significant) %>% 
      select(Gene_set)
    
    inter <- results.subset[[i]] %>% filter(ANOVA_group == "inter") %>% 
      mutate(Rank = min_rank(PValue)) %>% 
      filter(Rank <= n.significant) %>% 
      select(Gene_set)
    
    f <- results.subset[[i]] %>% filter(ANOVA_group == "F") %>% 
      mutate(Rank = min_rank(PValue)) %>% 
      filter(Rank <= n.significant) %>% 
      select(Gene_set)
    
    top.significant.sets <- rbind(treat, toler, inter, f)
    unique.sets <- append(unique.sets, top.significant.sets[["Gene_set"]]) %>% unique()
    
    plot.table <- rbind(plot.table, results.subset[[i]] %>% filter(Gene_set %in% unique.sets))
  }
  
  plot.table$Database <- factor(plot.table$Database, levels = c("KEGG", "BP", "MF", "CC"))
  plot.table$Gene_set <- sub(".....", "", plot.table$Gene_set)
  return(plot.table)
}

# Run the above function to create tables containing all rows to be incorporated into plot as dots.

# For each analysis function (cam, camPR, and gst), there will be 2 plot tables. 1 will contain only the significant 
# rows (note: if a gene set was in the top 4 significant sets for tolerant but for interaction it was significant and 
# not top 4, it will still be included here.)

# Another table will contain significant and nonsignificant rows (ex: if a gene set was in the top 4 significant sets 
# for tolerant but was not significant at all for interaction, it will still be included here.)

# The 1st plot table described above is a subset of the 2nd plot table (for when PValue <= 0.05)
cam.plot.table <- create_plot_table(results.list.reformatted, "cam", 4)
cam.plot.table.only.sig <- cam.plot.table %>% filter(PValue <= 0.05)

camPR.plot.table <- create_plot_table(results.list.reformatted, "camPR", 4)
camPR.plot.table.only.sig <- camPR.plot.table %>% filter(PValue <= 0.05)

gst.plot.table <- create_plot_table(results.list.reformatted, "gst", 4)
gst.plot.table.only.sig <- gst.plot.table %>% filter(PValue <= 0.05)

#--------------------------------------------------------------------------------------------------------------------------------
#Make dotplots again
  # For the cam/camPR dotplots, we will also incorporate the direction of significance using two separate color gradients.
create_camera_plots <- function(data.to.plot){
    # Sets the order for the x-axis in the plots. 
  x.axis <- factor(data.to.plot$ANOVA_group, levels = c("treat", "toler", "inter", "F"), 
                   labels = c("Treatment", "Tolerance", "Interaction", "F"))
  x.axis.down <- factor(data.to.plot[data.to.plot$Direction == "Down", "ANOVA_group"], 
                        levels = c("treat", "toler", "inter", "F"), 
                        labels = c("Treatment", "Tolerance", "Interaction", "F"))
    # Create the plot and return.
  plot <- ggplot() + 
    geom_point(data.to.plot, mapping = aes(x = x.axis, y = Gene_set, size = Matches, color = PValue)) + 
    facet_grid(Database ~ ., space = "free_y", scales = "free") +
    scale_color_gradientn("P-value (Up)", colors = c("green4", "gray80"), limits = c(0, 0.05)) +
    new_scale_color() +
    geom_point(data.to.plot %>% filter(Direction == "Down"), mapping = aes(x = x.axis.down, y = factor(Gene_set), 
                                                                             size = Matches, color = PValue)) + 
    facet_grid(Database ~ ., space = "free_y", scales = "free") +
    scale_color_gradientn("P-value (Down)", colors = c("red4", "gray80"), limits = c(0, 0.05)) +
    theme_bw() + 
    labs(x = "ANOVA effect", y = "Gene set", size = 'Number of genes') + 
    scale_y_discrete(label = function(Gene_set) str_trunc(Gene_set, 35)) +
    guides(size = guide_legend(order = 1))
  
  return(plot)
}

  # Run function to get both camera plots
cam.plot <- create_camera_plots(cam.plot.table)
cam.plot.only.sig <- create_camera_plots(cam.plot.table.only.sig)

camPR.plot <- create_camera_plots(camPR.plot.table)
camPR.plot.only.sig <- create_camera_plots(camPR.plot.table.only.sig)


  # Create gst dotplot
create_gst_plot <- function(data.to.plot){
  x.axis <- factor(data.to.plot$ANOVA_group, levels = c("treat", "toler", "inter", "F"), 
                   labels = c("Treatment", "Tolerance", "Interaction", "F"))
  
  plot <- ggplot(data = data.to.plot, aes(x = x.axis, y = Gene_set, size = Matches, color = PValue)) +
    facet_grid(Database ~ ., space = "free_y", scales = "free") +
    geom_point() + 
    scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + 
    theme_bw() +
    labs(x = "ANOVA effect", y = "Gene set", color = 'P-value', size = 'Number of genes') +
    scale_y_discrete(label = function(Gene_set) str_trunc(Gene_set, 25)) +
    guides(size = guide_legend(order = 1))
  
  return(plot)
}

# Run function to create plots
gst.plot <- create_gst_plot(gst.plot.table)
gst.plot.only.sig <- create_gst_plot(gst.plot.table.only.sig)

  # Save files
ggsave(plot = cam.plot, filename = "camera_dotplot.jpeg", device = "jpeg", path = output.directory)
ggsave(plot = cam.plot.only.sig, filename = "camera_dotplot_only_sig.jpeg", device = "jpeg", path = output.directory)
ggsave(plot = camPR.plot, filename = "cameraPR_dotplot.jpeg", device = "jpeg", path = output.directory)
ggsave(plot = camPR.plot.only.sig, filename = "cameraPR_dotplot_only_sig.jpeg", device = "jpeg", path = output.directory)
ggsave(plot = gst.plot, filename = "geneSetTest_dotplot.jpeg", device = "jpeg", path = output.directory)
ggsave(plot = gst.plot.only.sig, filename = "geneSetTest_dotplot_only_sig.jpeg", device = "jpeg", path = output.directory)


