#!/usr/bin/env Rscript

# R script to create search for significant interesting GO terms

#install.packages("eulerr")

# load librarys
library(stringr)

# turn off scientific notation
options(scipen = 999)

# set the working directory
workingDir <- "/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes"
setwd(workingDir)

# read in BP GO term enrichment results
treatmentGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/treatment_BP_sigGO_terms.csv")
toleranceGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/tolerance_BP_sigGO_terms.csv")
interactionGO <- read.csv("/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/DEAnalysis/Genotypes/GOAnalysis/interaction_BP_sigGO_terms.csv")

# import modules GO data
modulesNames <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$")
modulesNames <- str_remove(modulesNames, "_BP_sigGO_terms.csv.flt")
modulesFiles <- list.files(path="/Users/bamflappy/PfrenderLab/OLYM_dMelUV/KAP4/NCBI/GCF_021134715.1/Biostatistics/WGCNA/Genotypes/GOAnalysis_OLYM_30/", pattern = "_BP_sigGO_terms.csv.flt$", full.names = TRUE)
modulesGO <- lapply(modulesFiles, read.csv)
names(modulesGO) <- modulesNames

# combine GO data
combinedGO <- modulesGO
combinedGO[["Treatment"]] <- data.frame(GO.ID = treatmentGO$GO.ID)
combinedGO[["Tolerance"]] <- data.frame(GO.ID = toleranceGO$GO.ID)
combinedGO[["Interaction"]] <- data.frame(GO.ID = interactionGO$GO.ID)

# geneIDs with repair or other interesting GO terms
# The repair pathway GO terms we specifically explored in the analysis were DNA repair (GO:0006281), mismatch repair (MMR; GO:0006298), base excision repair (BER; GO:0006284), homologous recombination (HR; GO:0035825), nucleotide excision repair (NER; GO:0006289), intrastrand crosslink repair (ICL repair; GO:0036297), double strand break repair (DSBR; GO:0006302), single strand break repair (SSBR; GO:0000012).
repairTerms <- list(MMR = "GO:0006298",
                    BER = "GO:0006284",
                    HR = "GO:0035825",
                    NER = "GO:0006289",
                    ICLR = "GO:0036297",
                    DSBR = "GO:0006302",
                    SSBR = "GO:0000012",
                    PDR = "GO:0006290",
                    PR = "GO:0000719"
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

# The stress response terms included response to stress (GO:0006950), response to oxidative stress (GO:0006979), cellular response to oxidative stress (GO:0034599), 
# response to photooxidative stress (GO:0080183), detection of oxidative stress (GO:0070994), response to reactive oxygen species (GO:0000302), 
# cellular response to reactive oxygen (GO:0034614), response to hydroperoxide (GO:0033194).
stressTerms <- list(stress = "GO:0006950",
                    oxidative = "GO:0006979",
                    cellOxidative = "GO:0034599",
                    cellReactiveOxy = "GO:0034614",
                    regTransOxidative = "GO:0043556",
                    MAPKKK = "GO:1990315",
                    hydroperoxide = "GO:0071447",
                    senescence = "GO:0090403",
                    symbiont = "GO:0052164",
                    regCellOxidative = "GO:1900407",
                    photoOx = "GO:0080183",
                    detectOx = "GO:0070994",
                    ROS = "GO:0000302",
                    responseHyper = "GO:0033194"
)

# geneIDs with DNA damage response GO terms
# GO:0006974 DNA damage response
# GO:0008630 intrinsic apoptotic signaling pathway in response to DNA damage,  GO:0140861 DNA repair-dependent chromatin remodeling, GO:0141112 broken chromosome clustering
# GO:0009432 SOS response, GO:0006281 DNA repair, GO:0042770 signal transduction in response to DNA damage, GO:0043247 telomere maintenance in response to DNA damage
damageTerms <- list(DNADR = "GO:0006974",
                    intrinsic = "GO:0008630",
                    signal = "GO:0042770",
                    SOS = "GO:0009432",
                    telomere = "GO:0043247"
)

# geneIDs with DNA damage response GO terms
# GO:0006281 DNA repair
# GO:0036297 interstrand cross-link repair, GO:0006290 pyrimidine dimer repair, GO:0097196 Shu complex, GO:0009380 excinuclease repair complex
# GO:0006282 regulation of DNA repair, GO:0006284 base-excision repair, GO:0006289 nucleotide-excision repair, GO:0010213 non-photoreactive DNA repair
# GO:0006298 mismatch repair, GO:0006307 DNA dealkylation involved in DNA repair, GO:0045739 positive regulation of DNA repair
# GO:0000012 single strand break repair, GO:0000725 recombinational repair, GO:0046787 viral DNA repair, GO:0006302 double-strand break repair
# GO:0070914 UV-damage excision repair, GO:0000731 DNA synthesis involved in DNA repair, GO:0051103 DNA ligation involved in DNA repair
# GO:0043504 mitochondrial DNA repair, GO:0106300 protein-DNA covalent cross-linking repair, GO:0045004 DNA replication proofreading
# GO:0045738 negative regulation of DNA repair, GO:1990391 DNA repair complex, GO:0006301 postreplication repair
repairTerms_all <- list(DNAR = "GO:0006281",
                    ICLR = "GO:0036297",
                    PDR = "GO:0006290",
                    SC = "GO:0097196",
                    ERC = "GO:0009380",
                    RDNAR = "GO:0006282",
                    BER = "GO:0006284",
                    NER = "GO:0006289",
                    nonPhoto = "GO:0010213",
                    MMR = "GO:0006298",
                    dealkylation = "GO:0006307",
                    posRDNAR = "GO:0045739",
                    SSBR = "GO:0000012",
                    RR = "GO:0000725",
                    viralDNAR = "GO:0046787",
                    DSBR = "GO:0006302",
                    UVER = "GO:0070914",
                    synthesis = "GO:0000731",
                    ligation = "GO:0051103",
                    mito = "GO:0043504",
                    CLR = "GO:0106300",
                    proofreading = "GO:0045004",
                    negRDNAR = "GO:0045738",
                    DNARC = "GO:1990391",
                    postreplication = "GO:0006301"
)


# loop over each repair term
print("REPAIR")
# loop over each module
for(i in 1:length(combinedGO)){
  for(j in 1:length(repairTerms)){
    test <- combinedGO[[i]][["GO.ID"]]
    if (repairTerms[[j]] %in% test) {
      print(paste(names(combinedGO)[i], test[test == repairTerms[[j]]], names(repairTerms[j])), quote = FALSE)
    }
  }
}

# loop over each radiation term
print("RADIATION")
# loop over each module
for(i in 1:length(combinedGO)){
  for(j in 1:length(radiationTerms)){
    test <- combinedGO[[i]][["GO.ID"]]
    if (radiationTerms[[j]] %in% test) {
      print(paste(names(combinedGO)[i], test[test == radiationTerms[[j]]], names(radiationTerms[j])), quote = FALSE)
    }
  }
}

# loop over each stress term
print("STRESS")
# loop over each module
for(i in 1:length(combinedGO)){
  for(j in 1:length(stressTerms)){
    test <- combinedGO[[i]][["GO.ID"]]
    if (stressTerms[[j]] %in% test) {
      print(paste(names(combinedGO)[i], test[test == stressTerms[[j]]], names(stressTerms[j])), quote = FALSE)
    }
  }
}

# loop over each stress term
print("DAMAGE")
# loop over each module
for(i in 1:length(combinedGO)){
  for(j in 1:length(damageTerms)){
    test <- combinedGO[[i]][["GO.ID"]]
    if (damageTerms[[j]] %in% test) {
      print(paste(names(combinedGO)[i], test[test == damageTerms[[j]]], names(damageTerms[j])), quote = FALSE)
    }
  }
}

# loop over each stress term
print("ALL REPAIR")
# loop over each module
for(i in 1:length(combinedGO)){
  for(j in 1:length(repairTerms_all)){
    test <- combinedGO[[i]][["GO.ID"]]
    if (repairTerms_all[[j]] %in% test) {
      print(paste(names(combinedGO)[i], test[test == repairTerms_all[[j]]], names(repairTerms_all[j])), quote = FALSE)
    }
  }
}


