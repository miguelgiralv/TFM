library(SummarizedExperiment)
library(NetActivity)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(tibble)
library(dplyr)



path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr/"

complete_path=paste0(path,tissue_path)

#cogemos los tejidos que pueden tener relacion con el alzheimer
tissues_list<-list(
  "mashr_Artery_Aorta_predict.txt",                         
  "mashr_Artery_Coronary_predict.txt",                      
  "mashr_Artery_Tibial_predict.txt",                        
  "mashr_Brain_Amygdala_predict.txt",                       
  "mashr_Brain_Anterior_cingulate_cortex_BA24_predict.txt", 
  "mashr_Brain_Caudate_basal_ganglia_predict.txt",          
  "mashr_Brain_Cerebellar_Hemisphere_predict.txt",          
  "mashr_Brain_Cerebellum_predict.txt",                     
  "mashr_Brain_Cortex_predict.txt",                         
  "mashr_Brain_Frontal_Cortex_BA9_predict.txt",             
  "mashr_Brain_Hippocampus_predict.txt",                    
  "mashr_Brain_Hypothalamus_predict.txt",                   
  "mashr_Brain_Nucleus_accumbens_basal_ganglia_predict.txt",
  "mashr_Brain_Putamen_basal_ganglia_predict.txt",          
  "mashr_Brain_Substantia_nigra_predict.txt",               
  "mashr_Cells_EBV-transformed_lymphocytes_predict.txt",    
  "mashr_Heart_Atrial_Appendage_predict.txt",               
  "mashr_Heart_Left_Ventricle_predict.txt",                 
  "mashr_Whole_Blood_predict.txt"
)

files_tissue_list<-list()

for (i in 1:length(tissues_list)){
  path_tissue<-paste0(complete_path, tissues_list[[i]])
  files_tissue_list[[i]]<-path_tissue
}

# primero cargamos los df de predictxcan, que almacenaremos en una lista 
df_list_tissues <- list()
tissue_names <- list()

for (i in 1:length(files_tissue_list)) {
  tissue_df <- read.csv(files_tissue_list[[i]], header = TRUE, sep = "\t")
  tissue_names[[i]] <- sub(".*mashr_(.*)_predict\\.txt$", "\\1", files_tissue_list[[i]])
  df_list_tissues[[tissue_names[[i]]]] <- tissue_df
  print(paste(tissue_names[[i]], "procesado"))
}

save(df_list_tissues, file = paste0(path, "results/netactivity/df_list_tissues_curated.RData"))

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

save(score_tissue_list, file = paste0(path, "results/netactivity/se_list_curated.RData"))


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
  
  score_tissue_filtered <- SummarizedExperiment(
    assays = list(expression = filtered_assay_data),
    colData = colData(score_tissue),
    rowData = filtered_row_data
  )
  score_tissue_list[[tissue_names[[i]]]] <- score_tissue_filtered
  message(paste(tissue_names[[i]], "procesado, filas restantes:", nrow(filtered_assay_data)))
}

save(score_tissue_list, file = paste0(path, "results/netactivity/se_list3_curated.RData"))

# correccion en toda la tabla

all_pvalues <- data.frame()

for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[tissue_names[[i]]]]
  
  colData(score_tissue)$DiseaseStatus <- as.factor(colData(score_tissue)$DiseaseStatus)
  
  colData(score_tissue)$DiseaseStatus <- relevel(colData(score_tissue)$DiseaseStatus, ref = "control")
  
  tissue_name <- tissue_names[[i]]
  if (tissue_name %in% c("Testis","Prostate","Ovary","Uterus", "Vagina")) {
    model_score <- model.matrix(~ DiseaseStatus, colData(score_tissue))
  } else {
    model_score <- model.matrix(~ DiseaseStatus + Sex, colData(score_tissue))
  }
  
  
  fit <- lmFit(assay(score_tissue), model_score) %>% eBayes()
  
  topTab_tissue <- topTable(fit, coef = 2, n = Inf)
  
  topTab_tissue <- rownames_to_column(topTab_tissue, var = "GO")
  
  topTab_tissue$tissue<-tissue_names[[i]]
  
  all_pvalues <- rbind(all_pvalues, topTab_tissue)
  
  print(paste(tissue_names[[i]], "computado"))
}

all_pvalues$Bonferroni <- p.adjust(all_pvalues$P.Value, method = "bonferroni")  
all_pvalues$FDR <- p.adjust(all_pvalues$P.Value, method = "fdr")  

all_pvalues_Bonferroni<-all_pvalues %>% select(-"FDR")
all_pvalues_Bonferroni<-all_pvalues_Bonferroni[order(all_pvalues_Bonferroni$Bonferroni), ]

all_pvalues_FDR<-all_pvalues %>% select(-"Bonferroni")
all_pvalues_FDR <- all_pvalues_FDR[order(all_pvalues_FDR$FDR), ]

head(all_pvalues_FDR)

head(all_pvalues_Bonferroni)


write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni_curated.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR_curated.txt"), row.names = TRUE)


# correccion tejido por tejido (es la que da resultados)

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


all_pvalues_Bonferroni_2<-all_pvalues_2 %>% select(-"FDR")
all_pvalues_Bonferroni_2<-all_pvalues_Bonferroni_2[order(all_pvalues_Bonferroni_2$Bonferroni), ]

all_pvalues_FDR_2<-all_pvalues_2 %>% select(-"Bonferroni")
all_pvalues_FDR_2 <- all_pvalues_FDR_2[order(all_pvalues_FDR_2$FDR), ]

head(all_pvalues_FDR_2)
head(all_pvalues_Bonferroni_2)

# GO:0051481 Brain_Hippocampus
# GO:0070572 Heart_Left_Ventricle
# GO:0001516 Heart_Left_Ventricle

write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni_2.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR_2.txt"), row.names = TRUE)



## expresion de genes:
hippocampus<-score_tissue_list[["Brain_Hippocampus"]]
heart_left<-score_tissue_list[["Heart_Left_Ventricle"]]

data.frame(Expression = as.vector(assay(hippocampus["GO:0051481", ])),
           cases=hippocampus$DiseaseStatus) %>%
  ggplot(aes(x = cases, y = Expression, col = cases)) +
  geom_boxplot() +
  theme_bw() + 
  ylab("NetActivity scores") +
  ggtitle("Expresion de GO:0051481")

weights <- rowData(hippocampus)["GO:0051481", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")


data.frame(Expression = as.vector(assay(heart_left["GO:0070572", ])),
           cases=heart_left$DiseaseStatus) %>%
  ggplot(aes(x = cases, y = Expression, col = cases)) +
  geom_boxplot() +
  theme_bw() + 
  ylab("NetActivity scores") +
  ggtitle("Expresion de GO:0070572")

weights <- rowData(heart_left)["GO:0070572", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")


data.frame(Expression = as.vector(assay(heart_left["GO:0001516", ])),
           cases=heart_left$DiseaseStatus) %>%
  ggplot(aes(x = cases, y = Expression, col = cases)) +
  geom_boxplot() +
  theme_bw() + 
  ylab("NetActivity scores") +
  ggtitle("Expresion de GO:0001516")


weights <- rowData(heart_left)["GO:0001516", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")


#REGRESION LOGISTICA DE LOS GO


# GO:0051481 Brain_Hippocampus
# GO:0070572 Heart_Left_Ventricle
# GO:0001516 Heart_Left_Ventricle

tissue_name <- "Brain_Hippocampus"
go_term_of_interest <- "GO:0051481"

standardize_and_run_logistic_regression <- function(score_tissue_list, tissue_name, go_term_of_interest) {

  # extraemos el summarisedexperiment del tejido de interes
  score_tissue <- score_tissue_list[[tissue_name]]
  
  #transponemos su df
  assay_t <- as.data.frame(t(assay(score_tissue)))
  #cogemos los metadatos
  coldata <- as.data.frame(colData(score_tissue))
  
  # unimos el df con los metadatos 
  assay_t <- rownames_to_column(assay_t, "Sample")
  
  merged_data <- merge(coldata, assay_t, by = "Sample")
  
    if ("GO_term" %in% colnames(assay_t)) {
  #cogemos el geo term que especificamos
    merged_data <- merged_data[merged_data$GO_term == go_term_of_interest, ]
  }
  
  # cogemos las columnas con nuestro GO de interes 
  go_term_column <- colnames(assay_t)[which(colnames(assay_t) == go_term_of_interest)]
  
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
    Odds_Ratio = odds_ratios,
    CI_Lower = exp(conf_int[, 1]),
    CI_Upper = exp(conf_int[, 2])
  )
  
  # Return results
  return(list(
    c("este es el modelo"),
    summary=summary(model),
    c("este es el odds_ratio (95IC)"),
    odds_ratios = odds_df
  ))
}

# GO:0051481 Brain_Hippocampus
# GO:0070572 Heart_Left_Ventricle
# GO:0001516 Heart_Left_Ventricle


standardize_and_run_logistic_regression(score_tissue_list, "Heart_Left_Ventricle", "GO:0070572")
standardize_and_run_logistic_regression(score_tissue_list, "Heart_Left_Ventricle", "GO:0001516")
standardize_and_run_logistic_regression(score_tissue_list, "Brain_Hippocampus", "GO:0051481")

#cargamos el archivo de plink de scores de prs de PGS004146
# meter el sexo en el modelo tambien

PRS <- read.table(paste0(path, "results/PRS/PGS004146_hmPOS_GRCh37/output_PGS004146_5.sscore"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

head(PRS)

PRS$PHENO1 <- as.factor(PRS$PHENO) 
PRS$SCORE1_AVG <- as.numeric(PRS$SCORE)  

columns_to_remove <- c("X.FID", "IID", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM")

PRS <- PRS[ , !(colnames(PRS) %in% columns_to_remove)]
filtered_data$DiseaseStatus <- ifelse(merged_data$DiseaseStatus == "control", 0, 1)
# binarizamos
PRS$PHENO1<-ifelse(PRS$PHENO1 == 1, 0, 1)
# estandarizamos
PRS$SCORE1_AVG<-scale(PRS$SCORE1_AVG)
model <- glm(PHENO1 ~ ., data = PRS, family = binomial)

summary(model)
#no es estadísticamente significativo
#calculamos odds_Ratios:
exp(coef(model))

coef_summary <- summary(model)$coefficients
# Compute 95% confidence intervals
conf_int <- confint(model,level = 0.95)  
# Convert to data frames
odds_ratios <- exp(coef(model))  
odds_df <- data.frame(
  Odds_Ratio = odds_ratios,
  CI_Lower = exp(conf_int[, 1]),
  CI_Upper = exp(conf_int[, 2])
)


#cargamos el archivo de plink de scores de prs de PGS000054
# meter el sexo en el modelo tambien

PRS2 <- read.table(paste0(path, "results/PRS/PGS000054_hmPOS_GRCh37/output_PGS000054.sscore"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Specify the path to your file
path <- "/path/to/your/directory/"

# Read the file with the correct headers
PRS2 <- read.table(paste0(path, "results/PRS/PGS000054_hmPOS_GRCh37/output_PGS000054.sscore"), 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")


# Display the first few rows of the DataFrame
head(PRS2)


PRS2$PHENO1 <- as.factor(PRS2$PHENO) 
PRS2$SCORE1_AVG <- as.numeric(PRS2$SCORE)  

columns_to_remove <- c("X.FID", "IID", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM")

PRS2 <- PRS2[ , !(colnames(PRS2) %in% columns_to_remove)]

PRS2$SCORE1_AVG<-scale(PRS2$SCORE1_AVG)
model2 <- glm(PHENO1 ~ ., data = PRS2, family = binomial) 
# ¿METER TAMBIEN EL SEXO?
# Extract coefficients and their confidence intervals
summary(model2)
exp(coef(model2))
#no es estadísticamente significativo

coef_summary <- summary(model2)$coefficients
# Compute 95% confidence intervals
conf_int <- confint(model2,level = 0.95)  
# Convert to data frames
odds_ratios <- exp(coef(model2))  
odds_df <- data.frame(
  Odds_Ratio = odds_ratios,
  CI_Lower = exp(conf_int[, 1]),
  CI_Upper = exp(conf_int[, 2])
)






