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
library(DESeq2)
library(airway)
data(airway)

library(SummarizedExperiment)


## -----------------------------------------------------------------------------
ddsSE <- DESeqDataSet(airway, design = ~ cell + dex)
vst <- varianceStabilizingTransformation(ddsSE)

## -----------------------------------------------------------------------------
out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
out

## ----warning = FALSE, message = FALSE-----------------------------------------
library(tidyverse)

## ----plot, fig.cap = "Standardization. Effect of standardization on gene expression values. Before represent the gene expression values passed to prepareSummarizedExperiment. After are the values obtained after prepareSummarizedExperiment."----
rbind(assay(vst[1:5, ]) %>%
        data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        gather(Sample, Expression, 1:8) %>%
        mutate(Step = "Before"),
      assay(out[1:5, ]) %>%
        data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        gather(Sample, Expression, 1:8) %>%
        mutate(Step = "After")) %>%
  mutate(Step = factor(Step, levels = c("Before", "After"))) %>%
  ggplot(aes(x = Gene, y = Expression, col = Gene)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ Step, scales = "free") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())        

## -----------------------------------------------------------------------------
assay(out["ENSG00000278637", ]) 

## ----message = FALSE, warning = FALSE-----------------------------------------
library(limma)
library(Fletcher2013a)
data(Exp1)
Exp1

#cargamos la expresion en el tejido adiposo subcutaneo como un df en R:
df_adipose_en <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_en/en_Adipose_Subcutaneous_predict.txt", header = TRUE, sep = "\t")

head(df_adipose_en)

dim(df_adipose_en)

individuals<-df_adipose_en[,1]

genes <- df_adipose_en[, -c(1, 2)]

colnames(genes) <- sub("\\.\\d+$", "", colnames(genes))


row_data <- data.frame(Gene = colnames(genes))
col_data <- data.frame(SampleID = individuals)

se_adipose_en <- SummarizedExperiment(
  assays = list(expression = t(as.matrix(genes))),  
  rowData = row_data,
  colData = col_data
)

out_adipose_en <- prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg")



# Warning message:
#In prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg") :
# 5031 genes present in the model not found in input data. The expression of all samples will be set to 0.



df_adipose_mash <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_mashr/mashr_Adipose_Subcutaneous_predict.txt", header = TRUE, sep = "\t")

dim(df_adipose_mash)




se



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
library(DESeq2)
library(airway)
data(airway)

## -----------------------------------------------------------------------------
ddsSE <- DESeqDataSet(airway, design = ~ cell + dex)
vst <- varianceStabilizingTransformation(ddsSE)

## -----------------------------------------------------------------------------
out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
out

## ----warning = FALSE, message = FALSE-----------------------------------------
library(tidyverse)

## ----plot, fig.cap = "Standardization. Effect of standardization on gene expression values. Before represent the gene expression values passed to prepareSummarizedExperiment. After are the values obtained after prepareSummarizedExperiment."----
rbind(assay(vst[1:5, ]) %>%
        data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        gather(Sample, Expression, 1:8) %>%
        mutate(Step = "Before"),
      assay(out[1:5, ]) %>%
        data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        gather(Sample, Expression, 1:8) %>%
        mutate(Step = "After")) %>%
  mutate(Step = factor(Step, levels = c("Before", "After"))) %>%
  ggplot(aes(x = Gene, y = Expression, col = Gene)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ Step, scales = "free") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())        

## -----------------------------------------------------------------------------
assay(out["ENSG00000278637", ]) 

## ----message = FALSE, warning = FALSE-----------------------------------------
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
