library(SummarizedExperiment)
library(NetActivity)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(tibble)


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
  
  score_tissue_filtered <- SummarizedExperiment(
    assays = list(expression = filtered_assay_data),
    colData = colData(score_tissue),
    rowData = filtered_row_data
  )
  score_tissue_list[[tissue_names[[i]]]] <- score_tissue_filtered
  message(paste(tissue_names[[i]], "procesado, filas restantes:", nrow(filtered_assay_data)))
}

save(score_tissue_list, file = paste0(path, "results/netactivity/se_list3.RData"))

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


write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR.txt"), row.names = TRUE)


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


all_pvalues_Bonferroni_2<-all_pvalues_2 %>% select(-"FDR")
all_pvalues_Bonferroni_2<-all_pvalues_Bonferroni_2[order(all_pvalues_Bonferroni_2$Bonferroni), ]

all_pvalues_FDR_2<-all_pvalues_2 %>% select(-"Bonferroni")
all_pvalues_FDR_2 <- all_pvalues_FDR_2[order(all_pvalues_FDR_2$FDR), ]

head(all_pvalues_FDR_2)
head(all_pvalues_Bonferroni_2)


write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni_2.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR_2.txt"), row.names = TRUE)



## expresion de genes:


esophagus<-score_tissue_list[["Esophagus_Muscularis"]]
hippocampus<-score_tissue_list[["Brain_Hippocampus"]]
heart_left<-score_tissue_list[["Heart_Left_Ventricle"]]

data.frame(Expression = as.vector(assay(esophagus["GO:0002467", ])),
           cases=esophagus$DiseaseStatus) %>%
  ggplot(aes(x = cases, y = Expression, col = cases)) +
  geom_boxplot() +
  theme_bw() + 
  ylab("NetActivity scores") +
  ggtitle("Expresion de GO:0002467")

weights <- rowData(esophagus)["GO:0002467", ]$Weights_SYMBOL[[1]]
data.frame(weight = weights, gene = names(weights)) %>%
  mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
  ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Weight") +
  xlab("Gene")


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



#REGRESION LOGISTICA TEJIDO A TEJIDO

all_coef <- data.frame()


for (i in 1:length(score_tissue_list)) {
  score_tissue <- score_tissue_list[[i]]
  
  # transponemos los datos del tejido 
  assay_t <- as.data.frame(t(assay(score_tissue)))
  
  coldata <- as.data.frame(colData(score_tissue))
  
  # cada fila es una muestra
  assay_t <- rownames_to_column(assay_t, "Sample")
  
  # unimos con coldata que tiene los metadatos de las muestras
  merged_data <- merge(coldata, assay_t, by = "Sample")
  
  # quitamos sexo y muestra 
  merged_data <- merged_data[, !(names(merged_data) %in% c("Sample", "Sex"))]
  
  # binarizamos 
  merged_data$DiseaseStatus <- ifelse(merged_data$DiseaseStatus == "control", 0, 1)
  
  model <- glm(DiseaseStatus ~ ., data = merged_data, family = binomial)
  
  # Extraemos el oddsratio
  odds_ratios<- exp(coef(model))
  
  # lo convertimos en un df
  odds_df <- as.data.frame(odds_ratios)
  odds_df <- rownames_to_column(odds_df, "GO_term")
  names(odds_df)[2] <- "Odds_Ratio"
  
  # añadimos el nombre del tejido
  odds_df$tissue <- tissue_names[[i]]
  
  # aseguramos que Odds_Ratio es numerico para ordenarlo
  odds_df$Odds_Ratio <- as.numeric(odds_df$Odds_Ratio)
  
  
  odds_df <- odds_df[order(odds_df$Odds_Ratio, decreasing = FALSE), ]
  
  # lo unimos al df que almacena todos los odds ratio
  all_coef <- rbind(all_coef, odds_df)
  
  print(paste(tissue_names[[i]], "computado"))
}


colnames(all_coef) <- c("GO_term", "odds_ratio","tissue")

all_coef<-all_coef[order(all_coef$odds_ratio, decreasing = TRUE), ]


all_coef$GO_term<- as.numeric(odds_df$Odds_Ratio)

head(all_coef)

write.table(all_coef, file = paste0(path,"results/netactivity/", "odds_ratio.txt"), row.names = FALSE, sep = "\t",quote = FALSE)



#cargamos el archivo de plink de scores de prs

PRS <- read.table(paste0(path, "results/PRS/PGS004146_hmPOS_GRCh37.profile"), header = TRUE, sep = "", stringsAsFactors = FALSE)

head(PRS)

PRS$PHENO <- as.factor(PRS$PHENO) 
PRS$SCORE <- as.numeric(PRS$SCORE)  

PRS <- PRS %>% select(-c(FID, IID, CNT, CNT2))


model <- glm(PHENO ~ ., data = PRS, family = binomial)

summary(model)
#no es estadísticamente significativo
#calculamos odds_Ratios:
exp(coef(model))

