library(SummarizedExperiment)
library(NetActivity)
library(limma)
library(SummarizedExperiment)
library(tidyverse)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr"

complete_path=paste0(path,tissue_path)


files_tissue_list<- list.files(path = complete_path,
                               pattern = "^mashr_.*_predict\\.txt$",
                               full.names = TRUE)

# primero cargamos los df de predictxcan, que almacenaremos en una lista 
df_list_tissues <- list()
tissue_names <- list()

for (i in 1:length(files_tissue_list)) {
  tissue_df <- read.csv(files_tissue_list[[i]], header = TRUE, sep = "\t")
  tissue_names[[i]] <- sub(".*mashr_(.*)_predict\\.txt$", "\\1", files_tissue_list[[i]])
  df_list_tissues[[tissue_names[[i]]]] <- tissue_df
  print(paste(tissue_names[[i]], "procesado"))
}

save(df_list_tissues, file = paste0(path, "results/netactivity/df_list_tissues.RData"))

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
  
  se_tissue <- SummarizedExperiment(
    assays = list(expression = as.matrix(transposed_genes)), # se guarda como matriz
    colData = col_data,
    rowData = row_data
  )
  
  out_tissue <- prepareSummarizedExperiment(se_tissue, "gtex_gokegg")
  
  scores_tissue <- computeGeneSetScores(out_tissue, "gtex_gokegg")
  
  score_tissue_list[[tissue_names[[i]]]] <- scores_tissue  
  
  print(paste(tissue_names[[i]], "computado"))
}

save(score_tissue_list, file = paste0(path, "results/netactivity/se_list2.RData"))

load(paste0(path, "results/netactivity/se_list2.RData"))


# ahora eliminaremos los GO que tengan un valor 0 en todos los individuos (ya que significa que netactivity no ha tenido como input alguno de estos genes)
zero_sum_rows_list <- list()

for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[tissue_names[[i]]]]
  assay_data <- assay(score_tissue)
  row_sums <- rowSums(assay_data)
  
  zero_sum_rows <- which(row_sums == 0)
  zero_sum_row_names <- rownames(assay_data)[zero_sum_rows]
  zero_sum_rows_list[[tissue_names[[i]]]] <- zero_sum_row_names
  message(paste(tissue_names[[i]], "filas con valor de cero:", length(zero_sum_row_names)))
  filtered_assay_data <- assay_data[!rownames(assay_data) %in% zero_sum_row_names, ]
  filtered_row_data <- rowData(score_tissue)[!rownames(rowData(score_tissue)) %in% zero_sum_row_names, ]

  # Create a new SummarizedExperiment object with the filtered data
  score_tissue_filtered <- SummarizedExperiment(
    assays = list(expression = filtered_assay_data),
    colData = colData(score_tissue),
    rowData = filtered_row_data
  )
  # Update the list with the filtered SummarizedExperiment object
  score_tissue_list[[tissue_names[[i]]]] <- score_tissue_filtered
  # Print the number of remaining rows after filtering
  message(paste(tissue_names[[i]], "procesado, filas restantes:", nrow(filtered_assay_data)))
}

save(score_tissue_list, file = paste0(path, "results/netactivity/se_list3.RData"))


all_pvalues <- data.frame()
all_info <- data.frame()

for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[tissue_names[[i]]]]
  
  # Set "control" as the reference level
  colData(score_tissue)$DiseaseStatus <- as.factor(colData(score_tissue)$DiseaseStatus)
  
  colData(score_tissue)$DiseaseStatus <- relevel(colData(score_tissue)$DiseaseStatus, ref = "control")
  
  tissue_name <- tissue_names[[i]]
  if (tissue_name %in% c("Testis","Prostate","Ovary","Uterus", "Vagina")) {
    model_score <- model.matrix(~ DiseaseStatus, colData(score_tissue))
  } else {
    model_score <- model.matrix(~ DiseaseStatus + Sex, colData(score_tissue))
  }
  
  
  fit <- lmFit(assay(score_tissue), model_score) %>% eBayes()
  
  fit_p_val<-fit$p.value
  
  fit_p_val_df <- as.data.frame(fit_p_val)

  fit_p_val_df$tissue <- tissue_names[[i]]
  fit_p_val_df <- fit_p_val_df %>% select("DiseaseStatusAD", "tissue")
  

  all_pvalues <- rbind(all_pvalues, fit_p_val_df)
  
  print(paste(tissue_names[[i]], "computado"))
}

all_pvalues$Bonferroni <- p.adjust(all_pvalues$DiseaseStatusAD, method = "bonferroni")  
all_pvalues$FDR <- p.adjust(all_pvalues$DiseaseStatusAD, method = "fdr")  

all_pvalues_Bonferroni<-all_pvalues %>% select("DiseaseStatusAD", "tissue","Bonferroni")
all_pvalues_Bonferroni<-all_pvalues_Bonferroni[order(all_pvalues_Bonferroni$Bonferroni), ]

all_pvalues_FDR<-all_pvalues %>% select("DiseaseStatusAD", "tissue","FDR")
all_pvalues_FDR <- all_pvalues_FDR[order(all_pvalues_FDR$FDR), ]

write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR.txt"), row.names = TRUE)

