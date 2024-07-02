library(SummarizedExperiment)
library(NetActivity)
library(limma)
library(SummarizedExperiment)
library(tidyverse)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr"

complete_path=paste0(path,tissue_path)

individuals_cases<-

files_tissue_list<- list.files(path = complete_path,
                               pattern = "^mashr_.*_predict\\.txt$",
                               full.names = TRUE)


df_list_tissues<-list()
tissue_names<-list()

for (i in 1:length(files_tissue_list))
{
  tissue_df<-read.csv(files_tissue_list[[i]], header = TRUE, sep = "\t")
  tissue_names[[i]]<-sub(".*mashr_(.*)_predict\\.txt$", "\\1", files_tissue_list[[i]])
  df_list_tissues[[tissue_names[[i]]]]<-tissue_df
  print(paste(tissue_names[[i]],"procesado"))
}  


individuals_cases <- read.csv(paste0(path,"results/metadata_gsm/merged_metadata.txt"), header = TRUE, sep = "\t")

score_tissue_list <- list()

for (i in 1:length(df_list_tissues)) {
  # Load each data frame
  tissue_df <- df_list_tissues[[i]]
  
  # Extract and clean individual IDs
  individuals <- tissue_df[, 1]
  individuals <- sub("^GSE33528_", "", individuals)
  tissue_df[, 1] <- individuals
  
  # Remove version numbers from gene names
  colnames(tissue_df) <- sub("\\.\\d+$", "", colnames(tissue_df))
  
  # Merge with metadata
  merged_individuals <- merge(
    x = data.frame(Individual.ID = individuals), 
    y = individuals_cases[, c("Individual.ID", "Disease.status", "Sex")],
    by = "Individual.ID",
    all.x = TRUE
  )
  
  # Select gene columns and exclude individual id ("IID")
  tissue_df <- tissue_df[, !names(tissue_df) %in% c("IID"), drop = FALSE]
  
  transposed_genes <- t(tissue_df[,-1])
  
  colnames(transposed_genes) <- individuals
  
  # Add row_data and col_data
  row_data <- DataFrame(Gene = rownames(transposed_genes))
  
  col_data <- DataFrame(
    Sample = merged_individuals$Individual.ID,
    DiseaseStatus = merged_individuals$Disease.status,
    Sex = merged_individuals$Sex
  )
  
  # Create the SummarizedExperiment object
  se_tissue <- SummarizedExperiment(
    assays = list(expression = as.matrix(transposed_genes)), # save as matrix
    colData = col_data,
    rowData = row_data
  )
  
  out_tissue <- prepareSummarizedExperiment(se_tissue, "gtex_gokegg")
  
  scores_tissue <- computeGeneSetScores(out_tissue, "gtex_gokegg")
  
  score_tissue_list[[tissue_names[[i]]]] <- scores_tissue  
  
  print(paste(tissue_names[[i]], "computado"))
}



fit_tissue_list <- list()

bonferroni_dfs<-data.frame()
fdr_dfs<-data.frame()

for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[tissue_names[[i]]]]
  
  # Construct the model matrix
  model_score <- model.matrix(~ DiseaseStatus + Sex, colData(score_tissue))
  
  # Fit the linear model and apply empirical Bayes moderation
  fit <- lmFit(assay(score_tissue), model_score) %>% eBayes()
  

  topTab_bonferroni <- topTable(fit, coef = 1:2, n = Inf, adjust.method = "bonferroni")
  
    
  topTab_FDR <- topTable(fit, coef = 1:2, n = Inf, adjust.method = "fdr")

  
  topTab_bonferroni$tissue<-tissue_names[[i]]
  
  topTab_FDR$tissue<-tissue_names[[i]]
  
  bonferroni_dfs <- rbind(bonferroni_dfs, topTab_bonferroni)
  
  fdr_dfs <- rbind(fdr_dfs, topTab_FDR)
  
  
  # Store the results in the list
  fit_tissue_list[[tissue_names[[i]]]] <- list(
    fit = fit,
    topTab_bonferroni = topTab_bonferroni,
    topTab_FDR = topTab_FDR
  )
  
  print(paste(tissue_names[[i]], "computado"))
}

write.csv(bonferroni_dfs, file = paste0(path,"results/netactivity/", "bonferroni.txt"), row.names = TRUE)
write.csv(fdr_dfs, file = paste0(path,"results/netactivity/","FDR.txt"), row.names = TRUE)











