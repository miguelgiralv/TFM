library(ggplot2)
library(SummarizedExperiment)
library(tidyverse)
library(patchwork)
library(limma)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr/"

load(paste0(path,  "results/netactivity/linear_reg.RData"))

load(paste0(path,  "results/netactivity/se_list_curated_TODO.RData")) # netactivity ejecutar todo

# Comprobamos ahora cuales son los genes detectados en estos GOs y si entre ellos 
# están los que identificamos en el analisis diferencial de expresion
all_genes_list<-list()

for (i in  1:nrow(lin_analysis))
{
  GO<-lin_analysis$GO[i]
  tissue<-lin_analysis$tissue[i]
  tissue_score<-score_tissue_list[[tissue]]
  index<-which(rownames(tissue_score) == GO)
  gene_list<-rowData(tissue_score)[[4]][[index]]
  all_genes_list[[GO]] <- gene_list
}


all_genes_merge <- unique(unlist(lapply(all_genes_list, names)))


unique(all_genes_merge %in% c("MFAP3", "FAM114A2","CORO2A")) # los genes que detectamos en el analisis diferencial de expresion
# vemos que no están en ninguno de los GOs

# Comprobamos si en estos GOs interfieren genes con asociacione con alzheimer:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9966419/

previous_genes <- c(
  "APOE", "CD33", "CR1", "CLU", "PICALM", "EXOC3L2", "BIN1", "TOMM40", "ABCA7", 
  "MS4A6", "EPHA1", "CD2AP", "CUGBP2", "PVRL2", "DGKB", "GWA-10q23.1 (PCDH21, LRIT1, RGR)", 
  "HPCAL1", "GWA-9q21.33", "PPP1R3B", "PLD3", "TREM2", "APOC1", "CAMK1D", "FBXL13", 
  "PLXNA4", "FBXL7", "FRMD4A", "CELF1", "FERMT2", "SLC24A4-RIN3", "CNTNAP2", 
  "GWA-18p11.32", "GWA-12q24.23", "OSBPL6", "PTPRG", "PDCL3", "COBL", "SLC10A2", 
  "MLH3", "FNBP4", "CEACAM19", "CLPTM1", "ADAM10", "BCKDK", "KAT8", "VKORC1", "ACE",
  "TREML2", "PLCG2", "IL-34", "ADAMTS4", "INPP5D", "HESX1", "CLNK", "HS3ST1", 
  "HLA-DRB1", "ZCWPW1", "SPDYE3", "CNTNAP2", "CLU","PTK2B", "ECHDC3", "APH1B", 
  "BZRAP1-AS1", "SUZ12P1", "ALPK2", "AC074212.3", "CEACAM16", "BCL3", "MIR8085", 
  "CBLC", "BCAM", "PVRL2", "APOC1", "APOC1P1", "APOC4", "APOC2", "LINC00158", 
  "MIR155HG", "MIR155", "LINC00515", "MRPL39", "JAM2", "C2orf74", "ATG10", "MS4A6A", 
  "ABCB9", "ZNF815", "TRA2A", "MED30", "LPXN", "IRAK3", "N4BP2L2", "UQCC", 
  "APOBEC3F", "SFNa", "ESPN", "GNAI3", "C9orf72", "MTMR3", "NYAP1", "SPDYE3", 
  "ECHDC3", "USP6NL", "IQCK", "WWOX", "TRANK1", "FABP2", "LARP1B", "TSRM", "ARAP1", 
  "STARD10", "SPHK1", "SERPINB13", "EDEM1", "ALCAM", "GPC6", "VRK3", "SIPA1L2", 
  "WDR70", "API5", "ACER3", "PIK3C2G", "ARRDC4", "IGF1R", "RBFOX1", "MSX2", "AKAP9", 
  "PILRA", "NCK2", "SPI1", "TSPAN14", "SPPL2A", "ACE", "CCDC6", "ADAMTS1", "SHARPIN", 
  "GRN", "SPRED2", "ADAMTS4", "TMEM163", "SIGLEC11", "PLCG2", "IGHG1", "IKZF1", 
  "TSPOAP1", "FAM126A", "ZFHX4", "LGR5", "ZFC3H1", "OR51G1", "OR4X2", "ARHGEF4", 
  "PRUNE2", "MLKL", "NCOR2", "DMD", "NEDD4", "PLEC", "UNC5CL", "EPDR1", "WNT3", 
  "SORT1", "ADAM17", "PRKD3", "MME", "IDUA", "RHOH", "ANKH", "COX7C", "TNIP1", 
  "RASGEF1C", "HS3STS", "UMAD1", "ICA1", "TMEM106B", "JAZF1", "SEC61G", "CTSB", 
  "ABCA1", "ANK3", "BLNK", "PLEKHA1", "TPCN1", "SNX1", "DOC2A", "MAF", "FOXF1", 
  "PRDM7", "WDR91", "MYO15A", "KLF16", "LILRB2", "RBCK1", "SLC2A4RG", "APP"
)

additional_genes <- c(
  "ADAM10", "ABCA7", "BIN1", "CD2AP", "EPHA1", "PICALM", "SORL1", "MS4A6A", "CR1", 
  "CLU", "CD33", "TREM2", "TOMM40", "MAPT", "FERMT2", "CASS4", "PTK2B", "INPP5D", 
  "SLC10A2", "COBL", "UNC5C", "PLD3", "SLC24A4/RIN3", "HLA-DRB5","DRB1", "DSG2", 
  "MTHFR", "CST3", "BCHE", "CTSD", "ZCWPW1", "MEF2C", "ABI3", "PLCG2", "SCIMP", 
  "SHARPIN", "MINK1", "APH1B", "HS3ST1", "ECHDC3", "ACE", "PILRA", "SPI1", "IGF1", 
  "INSR", "LKB1"
)



combined_genes_bibliography <- unique(c(previous_genes, additional_genes))

common_genes <- intersect(all_genes_merge, combined_genes_bibliography)
common_genes
# "BCL3" "MEF2C" "ADAM17" "GRN" "IGF1R" 

# numero de genes totales que tomo como input netactivity:
all_genes_list<-list()

for (i in  1:nrow(lin_analysis))
{
  GO<-lin_analysis$GO[i]
  tissue<-lin_analysis$tissue[i]
  tissue_score<-score_tissue_list[[tissue]]
  index<-which(rownames(tissue_score) == GO)
  gene_list<-rowData(tissue_score)[[4]][[index]]
  all_genes_list[[GO]] <- gene_list
}


