---
  title: "Analysis of Subcutaneous Adipose Tissue Gene Expression"
author: "Your Name"
date: "`r format(Sys.time(), '%Y-%m-%d')`"
output:
  html_document:
  toc: true
toc_depth: 3

---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

## ----eval = FALSE-------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NetActivity")

## ----eval = FALSE-------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("yocra3/NetActivity")
#  devtools::install_github("yocra3/NetActivityData")

## ----warning = FALSE, message = FALSE-----------------------------------------
library(NetActivity)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(SparseArray)



#cargamos la expresion en el tejido adiposo subcutaneo como un df en R:
df_adipose_en <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_en/processed/processed_en_Adipose_Subcutaneous_predict.txt", header = TRUE, sep = "\t")

dim(df_adipose_en)


#procesamos las filas para que en los individuos no aparezca el gse:


individuals <- df_adipose_en[, 1]
individuals <- sub("^GSE33528_", "", individuals)

df_adipose_en[, 1] <- individuals



# tambien cambiamos los nombres de los genes para que no aparezca la version de ensembl
colnames(df_adipose_en) <- sub("\\.\\d+$", "", colnames(df_adipose_en))

# cargamos tambien los metadatos para saber que individuo tiene alzheimer y cual no
individuals_cases <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", header = TRUE, sep = "\t")
colnames(individuals_cases)

# fusionamos las dos matrices de individuos con metadatos y creamos una nueva que contenga solo los datos del experimento y su enfermedad
merged_individuals <- merge(
  x = data.frame(Individual.ID = individuals), 
  y = individuals_cases[, c("Individual.ID", "Disease.status")],
  by.x = "Individual.ID",
  by.y = "Individual.ID",
  all.x = TRUE  
)

# ahora transponemos nuestra matriz de expresion de forma que para cada fila haya un gen y para cada individuo haya una columna

gene_columns <- df_adipose_en[, !names(df_adipose_en) %in% "IID", drop = FALSE]

transposed_genes <- t(as.matrix(gene_columns))

df_adipose_en_trans <- as.data.frame(transposed_genes)


colnames(df_adipose_en_trans)<-individuals


df_adipose_en_trans[1,]

dim(df_adipose_en_trans)

#finalmente añadimos row_data y col_data

genes <- df_adipose_en[, !names(df_adipose_en) %in% "IID"]


row_data <- DataFrame(
  Gene = rownames(df_adipose_en_trans)
)

col_data <- data.frame(
  Sample = merged_individuals$Individual.ID,
  DiseaseStatus = merged_individuals$Disease.status
)



se_adipose_en <- SummarizedExperiment(
  assays = list(expression = as.matrix(df_adipose_en_trans)), # lo guardamos como una matriz y no un df
  colData = col_data,
  rowData = row_data
)
str(se_adipose_en)

# preparamos el experimento
out_adipose_en <- prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg")

# 5031 genes del modelo de netactivity no estan en nuestros datos, es decir que solo coge unos 3000 genes
scores_adipose_en <- computeGeneSetScores(out_adipose_en, "gtex_gokegg")
rownames(assay(scores_adipose_en))="adipose"

class(assay((scores_adipose_en)))

scores_adipose_en



assay(scores_adipose_en)
# le añade a rowData los pesos de los g
rowData(scores_adipose_en)
colData(scores_adipose_en)

mod_adipose_en <- model.matrix(~ DiseaseStatus , colData(scores_adipose_en))
fit <- lmFit(assay(scores_adipose_en), mod_adipose_en) %>% eBayes()
#Warning message:
#  Zero sample variances detected, have been offset away from zero
topTab <- topTable(fit, coef = 1:2, n = Inf)
head(topTab)




## -----------------------------------------------------------------------------
sessionInfo()
