library(SummarizedExperiment)
library(NetActivity)
library(limma)
library(tidyverse)
library(tibble)


path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr"

complete_path=paste0(path,tissue_path)


files_tissue_list<- list.files(path = complete_path,
                               pattern = "^mashr_.*_predict\\.txt$",
                               full.names = TRUE)

# primero cargamos los df de predixcan, que almacenaremos en una lista 
df_list_tissues <- list()
tissue_names <- list()

for (i in 1:length(files_tissue_list)) {
  tissue_df <- read.csv(files_tissue_list[[i]], header = TRUE, sep = "\t")
  tissue_names[[i]] <- sub(".*mashr_(.*)_predict\\.txt$", "\\1", files_tissue_list[[i]])
  df_list_tissues[[tissue_names[[i]]]] <- tissue_df
  print(paste(tissue_names[[i]], "procesado"))
}


# cargamos ahora nuestro archivo de metadatos
individuals_cases <- read.csv(paste0(path, "results/metadata_gsm/merged_metadata.txt"), header = TRUE, sep = "\t")

# y vamos a procesar nuestros dfs para convertirlos al objeto Summarisedexperiment para computarlos con Netactivity
score_tissue_list <- list()

for (i in 1:length(df_list_tissues)) {
  tissue_df <- df_list_tissues[[i]]
  
  # cambiamos la lista de individuos
  individuals <- tissue_df[, 1]
  individuals <- sub("^GSE33528_", "", individuals)
  tissue_df[, 1] <- individuals
  
  # cambiamos el nombre de los genes para que no aparezca la version
  colnames(tissue_df) <- sub("\\.\\d+$", "", colnames(tissue_df))
  
  # quitamos la columna IID ya que tenemos dos columnas de identificacio de individuos (FIID e IID)
  tissue_df <- tissue_df[, !names(tissue_df) %in% c("IID"), drop = FALSE]
  
  # creamos un unico df que incluye a los individuos de cada tejido y los metadatos, para que se tenga el mismo orden de individuos
  merged_individuals <- merge(
    x = data.frame(Individual.ID = individuals), 
    y = individuals_cases[, c("Individual.ID", "Disease.status", "Sex")],
    by = "Individual.ID",
    all.x = TRUE
  )
  # quitamos los individuos de tejidos que son especificos a un sexo 
  tissue_name <- tissue_names[[i]]
  if (tissue_name %in% c("Testis","Prostate")) {
    merged_individuals <- subset(merged_individuals, Sex == "male")
  } else if (tissue_name %in% c("Ovary","Uterus", "Vagina")) {
    merged_individuals <- subset(merged_individuals, Sex == "female")
  }
  
  tissue_df <- tissue_df[tissue_df[, 1] %in% merged_individuals$Individual.ID, , drop = FALSE]
  
  # ahora transformamos el df para que los genes sean filas en vez de columnas y las columnas sean los individuos
  transposed_genes <- t(tissue_df[,-1])
  colnames(transposed_genes) <- tissue_df[,1]
  
  # Finalmente creamos el objetoS ummarizedExperiment
  row_data <- DataFrame(Gene = rownames(transposed_genes))
  
  col_data <- DataFrame(
    Sample = merged_individuals$Individual.ID,
    DiseaseStatus = merged_individuals$Disease.status,
    Sex = merged_individuals$Sex
  )
  
  score_tissue <- SummarizedExperiment(
    assays = list(expression = as.matrix(transposed_genes)), # se guarda como matriz
    colData = col_data,
    rowData = row_data
  )
  score_tissue_list[[tissue_names[[i]]]] <- score_tissue  
  print(paste(tissue_names[[i]], "computado"))
}

save(score_tissue_list, file = paste0(path,  "results/genexp/Gene_expression.RData"))

load(paste0(path,  "results/genexp/Gene_expression.RData"))

# correccion tejido por tejido

all_pvalues_2 <- data.frame()

for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[tissue_names[[i]]]]
  
  colData(score_tissue)$DiseaseStatus <- as.factor(colData(score_tissue)$DiseaseStatus)
  
  colData(score_tissue)$DiseaseStatus <- relevel(colData(score_tissue)$DiseaseStatus, ref = "control")
  
  tissue_name <- tissue_names[[i]]
  
  model_score <- model.matrix(~ DiseaseStatus, colData(score_tissue))
  
  fit <- lmFit(assay(score_tissue), model_score) %>% eBayes()
  
  topTab_tissue <- topTable(fit, coef = 2, n = Inf)
  
  topTab_tissue <- rownames_to_column(topTab_tissue, var = "GO")
  
  topTab_tissue$tissue<-tissue_names[[i]]
  
  topTab_tissue$Bonferroni <- p.adjust(topTab_tissue$P.Value, method = "bonferroni")  
  
  topTab_tissue$FDR <- p.adjust(topTab_tissue$P.Value, method = "fdr")  
  
  all_pvalues_2 <- rbind(all_pvalues_2, topTab_tissue)
  
  print(paste(tissue_names[[i]], "computado"))
}

all_pvalues_3<-subset(all_pvalues_2, all_pvalues_2$FDR<0.05)

all_pvalues_3<-all_pvalues_3[order(all_pvalues_3$GO),]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = all_pvalues_3$GO,
  mart = ensembl
)
# METEMOS LOS NOMBRE DE LOS GENES
all_pvalues_4 <- merge(all_pvalues_3, gene_info, by.x = "GO", by.y = "ensembl_gene_id", all.x = TRUE)

names(all_pvalues_4)[names(all_pvalues_4) == "hgnc_symbol"] <- "Genes"

unique(all_pvalues_4$Genes)

complete_odds_df<-data.frame()

for (k in 1:nrow(all_pvalues_4))
{
  tissue_name<-all_pvalues_4$tissue[[k]]
  GO_term<-all_pvalues_4$GO[[k]]
  Genes<-all_pvalues_4$Genes[[k]]
  tissue_score<-score_tissue_list[[tissue_name]]
  
  go_term_of_interest<-paste0(GO_term,"_", tissue_name)
  print(go_term_of_interest)
  #transponemos su df
  assay_t <- as.data.frame(t(assay(tissue_score)))
  #cogemos los metadatos
  coldata <- as.data.frame(colData(tissue_score))
  
  # unimos el df con los metadatos 
  assay_t <- rownames_to_column(assay_t, "Sample")
  
  merged_data <- merge(coldata, assay_t, by = "Sample")
  
  # cogemos las columnas con nuestro GO de interes 
  go_term_column <- colnames(assay_t)[which(colnames(assay_t) == GO_term)]
  
  # quitamos el resto de columnas que no interesen
  filtered_data <- merged_data[, c("DiseaseStatus", "Sex", go_term_column)]
  
  #binarizar
  filtered_data$DiseaseStatus <- ifelse(merged_data$DiseaseStatus == "control", 0, 1)
  
  # Convertimos la expresion a valor numerico para ordenar 
  filtered_data[[go_term_column]] <- as.numeric(filtered_data[[go_term_column]])
  
  # Estandarizamos la expresión
  filtered_data[[go_term_column]] <- scale(filtered_data[[go_term_column]])
  
  # hacemos la regresion
  model <- glm(DiseaseStatus ~ ., data = filtered_data, family = binomial)
  
  coef_summary <- summary(model)$coefficients
  # Compute 95% confidence intervals
  conf_int <- confint(model,level = 0.95)  
  # Convert to data frames
  odds_ratios <- exp(coef(model))  
  odds_df <- data.frame(
    terms = c("Intercept", "Sex", GO_term),
    GO = GO_term,
    tissue=tissue_name,
    c_GO_term = go_term_of_interest,
    Odds_Ratio = odds_ratios,
    CI_Lower = exp(conf_int[, 1]),
    CI_Upper = exp(conf_int[, 2]),
    p_value=  coef_summary[,4])
  rownames(odds_df)<-NULL
  
  complete_odds_df <- rbind(complete_odds_df, odds_df)
}

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_info <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = complete_odds_df$GO,
  mart = ensembl
)
reg_analysis <- merge(complete_odds_df, gene_info, by.x = "GO", by.y = "ensembl_gene_id", all.x = TRUE)

names(reg_analysis)[names(reg_analysis) == "hgnc_symbol"] <- "Genes"


reg_analysis <- reg_analysis[!reg_analysis$terms %in% c("Intercept", "Sex"), ]


reg_analysis$signif <- ifelse(
  reg_analysis$p_value > 0.1, NA, 
  ifelse(
    reg_analysis$p_value < 0.1 & reg_analysis$p_value > 0.05, "*", 
    ifelse(
      reg_analysis$p_value < 0.05 & reg_analysis$p_value > 0.01, "**", 
      ifelse(
        reg_analysis$p_value < 0.01, "***", NA
      )
    )
  )
)


graph1 <- ggplot(reg_analysis, aes(x = Odds_Ratio, y = reorder(paste(GO, tissue, sep = "_"), Odds_Ratio), fill = Genes)) +
  # Main bars
  geom_bar(stat = "identity", position = "dodge", aes(x = Odds_Ratio - 1)) +
  
  # Error bars for confidence intervals
  geom_errorbar(aes(xmin = CI_Lower-1, xmax = CI_Upper-1), width = 0.2, position = position_dodge(0.9)) +
  
  # Significance annotations
  labs(title = "Odds Ratios for gene expression levels",
       x = "Odds Ratio",
       y = NULL,
       fill = "Genes") +
  # Adjust x-axis
  scale_x_continuous(limits = c(-0.5,1), breaks = seq(-0.5, 1, by = 0.25), labels = seq(0.5, 2, by = 0.25)) +
  
  # Reverse y-axis for correct label ordering
  scale_y_discrete(
    limits = rev(unique(paste(reg_analysis$GO, reg_analysis$tissue, sep = "_"))),
    labels = NULL) +
  theme(
    axis.title.x = element_text(size = 10),  # Adjust x-axis label size
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 12),
  )


hjust <- ifelse(reg_analysis$Odds_Ratio < 1, 1.25, -0.25) #para los valores menores de 1 le ponemos un hjust hacia la izqda (1.5)

graph2<-graph1 + 
  geom_text(
    aes(x = Odds_Ratio - 1, label = signif),
    hjust = hjust, 
    size = 4, 
    na.rm = TRUE, 
    position = position_dodge(1)
  )
graph2
# finalmente no usamos imagen
ggsave(paste0(path, "figures/genexp/gene_exp_reglog.png"), plot = graph2, width = 10, height = 3.5, dpi = 300)

