## ----setup, include=FALSE-----------------------------------------------------
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


#cargamos la expresion en el tejido adiposo subcutaneo como un df en R:
df_adipose_en <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_en/processed/processed_en_Adipose_Subcutaneous_predict.txt", header = TRUE, sep = "\t")

head(df_adipose_en)

dim(df_adipose_en)

individuals<-df_adipose_en[,1]

#procesamos las filas para que en los individuos no aparezca el gse:

first_column <- df_adipose_en[, 1, drop = FALSE]
first_column[, 1] <- sub("^GSE33528_", "", first_column[, 1])
df_adipose_en[, 1] <- first_column[, 1]

# tambien cambiamos los nombres de los genes para que no aparezca la version
colnames(df_adipose_en) <- sub("\\.\\d+$", "", colnames(df_adipose_en))

# cargamos tambien los metadatos para saber que individuo tiene alzheimer y cual no
individuals_cases <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", header = TRUE, sep = "\t")
colnames(individuals_cases)

# fusionamos las dos matrices de individuos con metadatos y creamos una nueva que contenga solo los datos del experimento y su enfermedad
merged_individuals <- merge(
  x = data.frame(Individual.ID = individuals),  # Ensure x is a dataframe with Individual.ID column
  y = individuals_cases[, c("Individual.ID", "Disease.status")],
  by.x = "Individual.ID",
  by.y = "Individual.ID",
  all.x = TRUE  # Keep all rows from x (individuals)
)

# ahora transponemos nuestra matriz de expresion de forma que para cada fila haya un gen y para cada individuo haya una columna

gene_columns <- df_adipose_en[, !names(df_adipose_en) %in% "IID", drop = FALSE]

transposed_genes <- t(as.matrix(gene_columns))

df_adipose_en_trans <- as.data.frame(transposed_genes)

colnames(df_adipose_en_trans)<-individuals

# Print or use new_df as needed
print(new_df)


assays = list(expression = t(as.matrix(df_adipose_en)))

expression[,1]

head(assays)
se_adipose_en <- SummarizedExperiment(
  assays,
  colData = col_data,
  rowData = row_data
)
assays = list(expression = t(as.matrix(df_adipose_en)))

vst_en <- varianceStabilizingTransformation(se_adipose_en)


out_adipose_en <- prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg")

scores_adipose_en <- computeGeneSetScores(out_adipose_en, "gtex_gokegg")

scores_adipose_en

rowData(scores_adipose_en)
colData(scores_adipose_en)

colData(scores)


mod_adipose_en <- model.matrix(~ Disease.status, colData(scores_adipose_en))
fit <- lmFit(assay(scores_adipose_en), mod_adipose_en) %>% eBayes()
topTab <- topTable(fit, coef = 1:2, n = Inf)
head(topTab)


----------------------------
library(limma)
library(Fletcher2013a)
data(Exp1)
Exp1

## -----------------------------------------------------------------------------
SE_fletcher <- SummarizedExperiment(exprs(Exp1), colData = pData(Exp1), rowData = fData(Exp1))

## ----warning=FALSE, message = FALSE-------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)

## -----------------------------------------------------------------------------
rownames(SE_fletcher) <- rowData(SE_fletcher)$SYMBOL
SE_fletcher <- SE_fletcher[!is.na(rownames(SE_fletcher)), ]
SE_fletcher <- SE_fletcher[!duplicated(rownames(SE_fletcher)), ]

rownames(SE_fletcher) <- mapIds(org.Hs.eg.db,
                                keys = rownames(SE_fletcher),
                                column = 'ENSEMBL',
                                keytype = 'SYMBOL')
SE_fletcher <- SE_fletcher[!is.na(rownames(SE_fletcher)), ]
SE_fletcher <- SE_fletcher[!duplicated(rownames(SE_fletcher)), ]
SE_fletcher

## -----------------------------------------------------------------------------
out_array <- prepareSummarizedExperiment(SE_fletcher, "gtex_gokegg")
out_array

## -----------------------------------------------------------------------------
scores <- computeGeneSetScores(out_array, "gtex_gokegg")
scores

## -----------------------------------------------------------------------------
rowData(scores)
colData(scores)

## -----------------------------------------------------------------------------
mod <- model.matrix(~ Treatment + Time, colData(scores))
fit <- lmFit(assay(scores), mod) %>% eBayes()
topTab <- topTable(fit, coef = 2:4, n = Inf)
head(topTab)

## -----------------------------------------------------------------------------
topTab$GeneSetName <- rowData(scores)[rownames(topTab), "Term"]
head(topTab)

## ----plotScores, fig.cap = "GO:1990440 activity scores. GO:1990440 presented the most significant difference due to treatment."----
data.frame(Expression = as.vector(assay(scores["GO:1990440", ])),
           Treatment = scores$Treatment) %>%
  ggplot(aes(x = Treatment, y = Expression, col = Treatment)) +
  geom_boxplot() +
  theme_bw() +
  ylab("NetActivity scores")

## ----plotWeight, fig.cap = "Weights of GO:1990440 gene set. The figure represents the weights used for computing the GO:1990440 gene set score. Weights are in absolute value to enable an easier comparison of their magnitude. Positive weights are shown in blue and negative in red.", fig.wide = TRUE----
weights <- rowData(scores)["GO:1990440", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")

## -----------------------------------------------------------------------------
sessionInfo()



## -----------------------------------------------------------------------------
SE_fletcher <- SummarizedExperiment(exprs(Exp1), colData = pData(Exp1), rowData = fData(Exp1))

## ----warning=FALSE, message = FALSE-------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)

## -----------------------------------------------------------------------------
rownames(SE_fletcher) <- rowData(SE_fletcher)$SYMBOL
SE_fletcher <- SE_fletcher[!is.na(rownames(SE_fletcher)), ]
SE_fletcher <- SE_fletcher[!duplicated(rownames(SE_fletcher)), ]

rownames(SE_fletcher) <- mapIds(org.Hs.eg.db,
                                keys = rownames(SE_fletcher),
                                column = 'ENSEMBL',
                                keytype = 'SYMBOL')
SE_fletcher <- SE_fletcher[!is.na(rownames(SE_fletcher)), ]
SE_fletcher <- SE_fletcher[!duplicated(rownames(SE_fletcher)), ]
SE_fletcher

## -----------------------------------------------------------------------------
out_array <- prepareSummarizedExperiment(SE_fletcher, "gtex_gokegg")
out_array

## -----------------------------------------------------------------------------
scores <- computeGeneSetScores(out_array, "gtex_gokegg")
scores

## -----------------------------------------------------------------------------
rowData(scores)

## -----------------------------------------------------------------------------
mod <- model.matrix(~ Treatment + Time, colData(scores))
fit <- lmFit(assay(scores), mod) %>% eBayes()
topTab <- topTable(fit, coef = 2:4, n = Inf)
head(topTab)

## -----------------------------------------------------------------------------
topTab$GeneSetName <- rowData(scores)[rownames(topTab), "Term"]
head(topTab)

## ----plotScores, fig.cap = "GO:1990440 activity scores. GO:1990440 presented the most significant difference due to treatment."----
data.frame(Expression = as.vector(assay(scores["GO:1990440", ])),
           Treatment = scores$Treatment) %>%
  ggplot(aes(x = Treatment, y = Expression, col = Treatment)) +
  geom_boxplot() +
  theme_bw() +
  ylab("NetActivity scores")

## ----plotWeight, fig.cap = "Weights of GO:1990440 gene set. The figure represents the weights used for computing the GO:1990440 gene set score. Weights are in absolute value to enable an easier comparison of their magnitude. Positive weights are shown in blue and negative in red.", fig.wide = TRUE----
weights <- rowData(scores)["GO:1990440", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")

## -----------------------------------------------------------------------------
sessionInfo()
