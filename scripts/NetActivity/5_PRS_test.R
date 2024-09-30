
library(tidyverse)
library(tibble)
library(dplyr)
library(pROC)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" 


load(paste0(path,  "results/netactivity/se_list_curated_TODO.RData"))

process_PRS <- function(PRS) 
{
  score_tissue<-score_tissue_list[[1]]
  # cogemos cualquier tejido, lo que importa es seleccionar los individuos del conjunto test
  assay_tissue<- assay(score_tissue)
  # guardamos los individuos seleccionados en este fold:
  samples_test<-colnames(assay_tissue)
  
  PRS_processed <- PRS[PRS[["IID"]] %in% samples_test, ]
  
  PRS_processed$PHENO1 <- as.factor(PRS_processed$PHENO1) 
  PRS_processed$PRS <- as.numeric(PRS_processed$SCORE1_AVG)  
  
  columns_to_remove <- c("X.FID","SCORE1_AVG", "IID", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM")
  
  PRS_processed <- PRS_processed[ , !(colnames(PRS_processed) %in% columns_to_remove)]
  # binarizamos
  PRS_processed$PHENO1<-ifelse(PRS_processed$PHENO1 == 1, 0, 1)
  # estandarizamos
  PRS_processed$PRS<-scale(PRS_processed$PRS)
  model <- glm(PHENO1 ~ ., data = PRS_processed, family = binomial)
  
  predicted_probs <- predict(model, PRS_processed, type = "response")
  # Clasificamos predicciones con umbral de 0.5
  predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
  
  roc_obj <- roc(PRS_processed$PHENO1, predicted_probs)
  auc_value <- auc(roc_obj)
  # Calculamos matriz confusion
  conf_matrix <- confusionMatrix(factor(predicted_classes), factor(PRS_processed$PHENO1))
  
  accuracy <- conf_matrix$overall['Accuracy']
  precision <- conf_matrix$byClass['Pos Pred Value']
  recall <- conf_matrix$byClass['Sensitivity']
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  #no es estadísticamente significativo
  #calculamos odds_Ratios:
  exp(coef(model))
  
  coef_summary <- summary(model)$coefficients
  # 95% CI
  conf_int <- confint(model,level = 0.95)  

    odds_ratios <- exp(coef(model))  
  odds_df <- data.frame(
    Odds_Ratio = odds_ratios,
    CI_Lower = exp(conf_int[, 1]),
    CI_Upper = exp(conf_int[, 2])
  )
  odds_df<-cbind(odds_df, coef_summary)
  
  metrics_df <- data.frame(
    Metric = c("AUC", "Accuracy", "Precision", "Recall", "F1 Score"),
    Value = c(auc_value, accuracy, precision, recall, f1_score)
  )
  return(list(odds_df,metrics_df))
}

#cargamos el archivo de plink de scores de prs de PGS004146

# leemos los archivos de score para PGS004146 (poblacion general)
PRS1 <- read.table(paste0(path, "results/PRS/PGS004146_hmPOS_GRCh37/output_PGS004146.sscore"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
# y el de poblacion hispana 
PRS2 <- read.table(paste0(path, "results/PRS/PGS000054_hmPOS_GRCh37/output_PGS000054.sscore"), 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

# aplicamos la funcion:

PRS_results_general<-process_PRS(PRS1)
PRS_results_hispana<-process_PRS(PRS2)

save(PRS_results_general, file = paste0(path,  "results/PRS1.RData"))
save(PRS_results_hispana, file = paste0(path,  "results/PRS2.RData"))