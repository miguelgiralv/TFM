library(SummarizedExperiment)
library(caret)
library(pROC)
library(limma)
library(tidyverse)
library(tibble)
library(dplyr)
library(glmnet)
library(stats)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

Metascore<-function(score_tissue_list,tissue_names)
{
  set.seed(300)
  score_tissue_test<-score_tissue_list[[1]]
  col_data <- as.data.frame(colData(score_tissue_test))
  folds <- createFolds(col_data$Sample, k = 5)
  score_tissue_folds <- list()
  
  # iteramos por cada fold
  for (fold in 1:5) 
  {
    # hacemos una lista para cada fold
    score_tissue_folds[[paste0("fold_", fold)]] <- list(train = list(), test = list())
  }
  # iteramos para cada tejido
  for (i in 1:length(tissue_names)) {
    tissue_name <- tissue_names[[i]]
    score_tissue <- score_tissue_list[[tissue_name]]
    # y en cada fold
    for (fold in 1:5) {
      sample_indices <- folds[[fold]]
      n_rows <- ncol(score_tissue)
      test_score <- score_tissue[, sample_indices]
      train_score <- score_tissue[, -sample_indices]
      score_tissue_folds[[paste0("fold_", fold)]][["train"]][[tissue_name]] <- train_score
      score_tissue_folds[[paste0("fold_", fold)]][["test"]][[tissue_name]] <- test_score
    }
  }
  # Ahora sacamos la regresion lineal en cada fold para el conjunto train
  all_pvalues_list <- list()
  for (fold in 1:5)
  {
    train_tissues <- score_tissue_folds[[fold]][["train"]]
    # guardaremos los pvalores de todos los tejidos en un df
    all_pvalues <- data.frame()
    # iteramos en todos los tejidos por cada fold
    for (tissue_name in names(train_tissues)) 
    {
      score_tissue <- train_tissues[[tissue_name]]
      colData(score_tissue)$DiseaseStatus <- as.factor(colData(score_tissue)$DiseaseStatus)
      colData(score_tissue)$DiseaseStatus <- relevel(colData(score_tissue)$DiseaseStatus, ref = "control")
      
      model_score <- model.matrix(~ DiseaseStatus, colData(score_tissue))
            
      fit <- lmFit(assay(score_tissue), model_score) %>% eBayes()
      
      topTab_tissue <- topTable(fit, coef = 2, n = Inf)
      topTab_tissue <- rownames_to_column(topTab_tissue, var = "GO")
      topTab_tissue$tissue <- tissue_name
      topTab_tissue$FDR <- p.adjust(topTab_tissue$P.Value, method = "fdr")
      
      # Añadimos los resultados 
      all_pvalues <- rbind(all_pvalues, topTab_tissue)
      print(paste(tissue_name, "fold ", fold," computado "))
    }
    # unimos lso pvalores de folds en l lista
    all_pvalues_list[[paste0("fold_", fold)]] <- all_pvalues
  }
  # ahora filtramos nuestras matrices de regresion lineal corregidas para cada fold
  # y añadimos una columna que une GO y tejido
  all_pvalues_list_filtered<-list()
  
  for (i in 1:length(all_pvalues_list))
  {
    fold_linearmodel<-all_pvalues_list[[i]]
    fold_linearmodel <- fold_linearmodel[fold_linearmodel$FDR < 0.25 &
                                           (fold_linearmodel$logFC > 0.05 | fold_linearmodel$logFC < -0.05), ]
    # eliminamos los GOs duplicados
    duplicate_rows <- fold_linearmodel %>%
      group_by(GO) %>%
      filter(n() > 1)
    
    filtered_unique <- fold_linearmodel %>%
      group_by(GO) %>%
      filter(abs(logFC) == max(abs(logFC))) %>%
      distinct(GO, .keep_all = TRUE) 
    
    filtered_unique<-as.data.frame(filtered_unique)
    filtered_unique$complete_GO<- paste0(filtered_unique$GO, "_", filtered_unique$tissue)
    all_pvalues_list_filtered[[paste0("fold_", i)]]<-filtered_unique
  }   
  ########### ahora hacemos la regresion logistica con los GOs seleccionados en cada fold para hacer el metascore
  # almacenaremos los gos y su tejido en orden en una lista
  
  merged_data_list <- list()

  for (fold in 1:length(all_pvalues_list_filtered)) 
  {
    filtered_go <- all_pvalues_list_filtered[[fold]]
    data_list <- list()
    
    for (i in 1:nrow(filtered_go)) 
    {
      GO_name <- filtered_go$GO[i]
      go_tissue_name <- filtered_go$tissue[i]
      complete_GO <- filtered_go$complete_GO[i]
      
      # Extraemos el summarisedexperiment para cada uno de los los tejidos de cada GO 
      score_tissue <- score_tissue_folds[[fold]][["train"]][[go_tissue_name]]
      
      # Transponemos lo tejidos y cogemos los GSAS para cada GO en su tejido especifico
      assay_t <- as.data.frame(t(assay(score_tissue)))
      data_list[[complete_GO]] <- assay_t[[GO_name]]
      print(paste("Fold", fold, "- Processed GO:", complete_GO))
    }
    
    # Combinamos en un unicoo df
    combined <- as.data.frame(data_list)
    colnames(combined)<-gsub("^GO\\.", "GO:", colnames(combined))
    colnames(combined) <- gsub("^path.", "path:", colnames(combined))
    colnames(combined)  <- gsub('\\.', '-', colnames(combined))
    coldata <- as.data.frame(colData(score_tissue))
    rownames(coldata) <- NULL
    
    # Combinamos los datos del assay con los metadatos
    merged_data <- cbind(combined, coldata)
    merged_data <- dplyr::select(merged_data,-"Sample")
    
    # Binarizamos la enfermedad 
    merged_data$DiseaseStatus <- ifelse(merged_data$DiseaseStatus == "control", 0, 1)
    
    exclude_cols <- c("DiseaseStatus", "Sex")
    for (col_name in colnames(merged_data)) 
    {
      if (!(col_name %in% exclude_cols)) 
        {
        merged_data[[col_name]] <- as.numeric((merged_data[[col_name]]))
        merged_data[[col_name]] <- scale(merged_data[[col_name]])
        }
    }
    # ahora hacemos la regresion logistica con cada fold
    model <- glm(DiseaseStatus ~ ., data = merged_data, family = binomial)
    coef_summary <- as.data.frame(summary(model)$coefficients)
    coef_summary<-rownames_to_column(coef_summary, "Predictors")
    coef_summary[["Pr(>|z|)"]]<-as.numeric(coef_summary[["Pr(>|z|)"]])                                 
    coef_summary[["Estimate"]]<-as.numeric(coef_summary[["Estimate"]])
    # ponemos un umbral de Pr(>|z|) de 0.10
    coef_summary <- coef_summary[coef_summary[["Pr(>|z|)"]] < 0.10 | coef_summary[["Predictors"]]=="(Intercept)", ] 

    # y lo almacenamos en la lista
    merged_data_list[[paste0("fold_", fold)]] <- coef_summary
    print(paste("Fold", fold, "completed and added to merged_data_list."))
  }
  
  ####AHORA HAREMOS EL METASCORE
  # primero extraer los scores de GSAS de los terminos GO para cada muestra
  Metascore_test_list <- list()
  
  for (fold in 1:5)
  {
    models_fold <- merged_data_list[[fold]]
    models_fold$Predictors <- gsub('\`', '', models_fold$Predictors)
    models_fold$Predictors <- gsub("^path.", "path:", models_fold$Predictors)
    models_fold$Predictors <- gsub('\\.', '-', models_fold$Predictors)
    fold_combine_score_list<-list()
    
    for (predictor in models_fold$Predictors[-1]) # nos saltamos el primer elemento que es el intercept
    {  
      # dividimos para sacar el GO y tejido
      split_result <- strsplit(predictor, "_", fixed = TRUE)[[1]]
      tissue_name <- paste(split_result[-1], collapse = "_")
      # corregimos posibles errores tipograficos acumulados
      GO_tissue_se<-score_tissue_folds[[fold]][["test"]][[tissue_name]]
      colData(GO_tissue_se)$DiseaseStatus <- as.factor(colData(GO_tissue_se)$DiseaseStatus)
      colData(GO_tissue_se)$DiseaseStatus <- relevel(colData(GO_tissue_se)$DiseaseStatus, ref = "control")
      assay_tissue<- assay(GO_tissue_se)
      rownames(assay_tissue) <- paste0(rownames(assay_tissue),"_", tissue_name)
      assay_tissue<-t(assay_tissue)
      assay_tissue<-as.data.frame(assay_tissue)
      rownames(assay_tissue) <- NULL
      fold_combine_score_list[[predictor]]<-assay_tissue[predictor]
    }
    fold_metascore_df<-as.data.frame(fold_combine_score_list)
    # ponemos los nombre bien
    colnames(fold_metascore_df) <- gsub('GO.', 'GO:', colnames(fold_metascore_df))
    
    # cogemos datos del ultimo bucle
    coldata <- as.data.frame(colData(GO_tissue_se))
    rownames(coldata) <- NULL
    
    # combinamos nuestro df con los metadatos quitando sexo y control
    fold_metascore_df <- cbind(fold_metascore_df, coldata)
    fold_metascore_df <- dplyr::select(fold_metascore_df,-"Sample",-"Sex")
    
    # Binarizamos
    fold_metascore_df$DiseaseStatus <- ifelse(fold_metascore_df$DiseaseStatus == "control", 0, 1)
  
    # escalamos con scale 
    exclude_cols <- c("DiseaseStatus")
    for (col_name in colnames(fold_metascore_df)) 
    {
      if (col_name %in% models_fold$Predictors)
      {
        #multiplicamos los de GSAS para cada GO en cada individuo por el valor continuo del estimate para ese GO procedente
        # de nuestra regresion logistica con GOs seleccionados (models_fold)
        estimate_value <- models_fold[models_fold[["Predictors"]] == col_name, "Estimate"]
        fold_metascore_df[[col_name]]<-fold_metascore_df[[col_name]]*estimate_value
      }
    }
    # añadimos el valor del intercept al df:
    intercept_value <- models_fold[models_fold[["Predictors"]] == "(Intercept)", "Estimate"]
    fold_metascore_df$intercept=-intercept_value
    
    fold_computed_metascore<-data.frame()
    #  sumamos todas las columnas (suma ponderada) para obtener el metascore analogo al poligenic risk score
    metascore_score<- rowSums(fold_metascore_df[, !colnames(fold_metascore_df) %in% "DiseaseStatus"])
    fold_metascore_df <- dplyr::select(fold_metascore_df,"DiseaseStatus")
    fold_metascore_df$metascore<-metascore_score
    
    # por ultimo haremos la regresion logistica de este fold:
    # cogemos solo metascore y disease status
    fold_metascore_df$DiseaseStatus <- as.factor(fold_metascore_df$DiseaseStatus) 
    fold_metascore_df$metascore <- as.numeric(fold_metascore_df$metascore)  
    fold_metascore_df$metascore <- scale(fold_metascore_df$metascore)  
  
    model_metascore <- glm(DiseaseStatus ~ ., data = fold_metascore_df, family = binomial) 
    coef_summary <- summary(model_metascore)$coefficients
    # Intervalos de confianza del 95% 
    conf_int <- confint(model_metascore,level = 0.95)  
    # lo metemos en df
    odds_ratios <- exp(coef(model_metascore))  
    odds_df <- data.frame(
      Odds_Ratio = odds_ratios,
      CI_Lower = exp(conf_int[, 1]),
      CI_Upper = exp(conf_int[, 2]))
    odds_df<-cbind(odds_df, coef_summary)
    Metascore_test_list[[paste0("fold_",fold)]]<-list(reg_log_model = model_metascore,
                                                      odds_ratio=odds_df, 
                                                      GOs_metascore_regression = models_fold, 
                                                      metascore = fold_metascore_df)
  }
return(Metascore_test_list)
}

# CON TODOS  LOS  TEJIDOS:
load(paste0(path,  "results/netactivity/se_list_curated_TODO.RData"))

# quitamos los tejidos sexo especificos
excluded_tissues <- c("Testis", "Uterus", "Vagina","Prostate", "Ovary")

tissue_names <- setdiff(names(score_tissue_list), excluded_tissues)

Metascore_test_list_all<-Metascore(score_tissue_list,tissue_names)

# con tejidos seleccionados:
load(paste0(path,  "results/netactivity/se_list_curated.RData"))
tissue_names <- names(score_tissue_list)

Metascore_test_list_selected<-Metascore(score_tissue_list,tissue_names)

# obervamos ahora los resultados para cada fold:
Metascore_test_list_all[["fold_1"]]$odds_ratio
Metascore_test_list_all[["fold_2"]]$odds_ratio
Metascore_test_list_all[["fold_3"]]$odds_ratio
Metascore_test_list_all[["fold_4"]]$odds_ratio
Metascore_test_list_all[["fold_5"]]$odds_ratio


# calculamos ahora las performance metrics de los dos metascores
performance_metrics_calculator<-function(Metascore_test_list)
{
performance_metrics <- list()

for (fold in 1:5) {
  fold_results <- Metascore_test_list[[paste0("fold_", fold)]]
  metascore_df <- fold_results$metascore
  # calculamos auc y roc
  roc_obj <- roc(metascore_df$DiseaseStatus, metascore_df$metascore)
  auc_value <- auc(roc_obj)
  reg_log_model<-Metascore_test_list[[paste0("fold_", fold)]][["reg_log_model"]]
  predicted_probs <- predict(reg_log_model, fold_results$metascore, type = "response")
  predicted <- ifelse(predicted_probs > 0.5, 1, 0)
  # Calculamos matriz de confusion
  confusion <- confusionMatrix(factor(predicted), factor(metascore_df$DiseaseStatus))
  # Extraemos las metricas
  accuracy <- confusion$overall["Accuracy"]
  precision <- confusion$byClass["Pos Pred Value"]
  recall <- confusion$byClass["Sensitivity"]
  # Y las metemos en un df para cada fold dentro de una lista
  performance_metrics[[paste0("fold_", fold)]] <- data.frame(
    Fold = fold,
    AUC = auc_value,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall
  )
}

# combinamos en un unico df
performance_metrics_df <- do.call(rbind, performance_metrics)

# Calculamnos tb la media 
average_metrics <- performance_metrics_df %>%
  summarise(
    Fold = "Average",
    AUC = mean(AUC),
    Accuracy = mean(Accuracy),
    Precision = mean(Precision),
    Recall = mean(Recall)
  )

metrics_summary <- rbind(performance_metrics_df, average_metrics)

metrics_summary <- metrics_summary %>%
  mutate(across(c(AUC, Accuracy, Precision, Recall), as.numeric))

Metascore_test_list[["Performance"]]<-metrics_summary
return(Metascore_test_list)
}

Metascore_test_list_all<-performance_metrics_calculator(Metascore_test_list_all)
Metascore_test_list_selected<-performance_metrics_calculator(Metascore_test_list_selected)

save(Metascore_test_list_all, file = paste0(path, "results/netactivity/Metascore_TODO.RData"))
save(Metascore_test_list_selected, file = paste0(path, "results/netactivity/Metascore_selected.RData"))
