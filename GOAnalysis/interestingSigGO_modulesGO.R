#!/usr/bin/env Rscript

# R script to create search for significant interesting GO terms

#install.packages("eulerr")

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
setwd(workingDir)

# import modules GO data
modulesNames <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv$")
modulesNames <- str_remove(modulesNames, "_BP_sigGO_terms.csv")
modulesFiles <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv$", full.names = TRUE)
modulesGO <- lapply(modulesFiles, read.csv)
names(modulesGO) <- modulesNames

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


# loop over each repair term
print("REPAIR")
# loop over each module
for(i in 1:length(modulesGO)){
  for(j in 1:length(repairTerms)){
    test <- modulesGO[[i]][["GO.ID"]]
    if (repairTerms[[j]] %in% test) {
      print(names(modulesGO)[i])
      print(test[test == repairTerms[[j]]])
    }
  }
}

# loop over each radiation term
print("RADIATION")
# loop over each module
for(i in 1:length(modulesGO)){
  for(j in 1:length(radiationTerms)){
    test <- modulesGO[[i]][["GO.ID"]]
    if (radiationTerms[[j]] %in% test) {
      print(names(modulesGO)[i])
      print(test[test == radiationTerms[[j]]])
    }
  }
}

# loop over each stress term
print("STRESS")
# loop over each module
for(i in 1:length(modulesGO)){
  for(j in 1:length(stressTerms)){
    test <- modulesGO[[i]][["GO.ID"]]
    if (stressTerms[[j]] %in% test) {
      print(names(modulesGO)[i])
      print(test[test == stressTerms[[j]]])
    }
  }
}
