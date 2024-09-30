library(ggplot2)
library(SummarizedExperiment)
library(tidyverse)
library(patchwork)
library(limma)
library(cowplot)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio

tissue_path="results/predictXcan/test/vcf_1000G_hg37_mashr/"

complete_path=paste0(path,tissue_path)

# ya hemos hecho el muestreo, ahora haremos una regresión lineal con los individuos muestreados: 
pvalue_boxplot<-function(score_tissue_list)
{  
  tissue_names <- names(score_tissue_list)
  
  all_pvalues <- data.frame()
  
  for (i in 1:length(score_tissue_list)) 
  {
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
    
    topTab_tissue$FDR <- p.adjust(topTab_tissue$P.Value, method = "fdr")  
    
    all_pvalues <- rbind(all_pvalues, topTab_tissue)
    
    print(paste(tissue_names[[i]], "computado"))
  }
  
  #ahora los ordenamos
  all_pvalues_2<-all_pvalues[all_pvalues$FDR< 0.05, ]
  
  # vemos cuanto tienen un p vlaor <0.05
  sum(all_pvalues_2$P.Value<0.05)
  
  # y cuantos tienen un FDR<0.05
  sum(all_pvalues_2$FDR<0.05)
  
  # REGRESION LOGISTICA Y BARPLOT
  
  plot_list <- list()
  
  y_axis_limit <- 0.37  
  
  for (j in 1:nrow(all_pvalues_2))
  {
    tissue_name<-all_pvalues_2$tissue[[j]]
    GO_term<-all_pvalues_2$GO[[j]]
    
    tissue_score<-score_tissue_list[[tissue_name]]
    weights <- rowData(tissue_score)[GO_term, ]$Weights_SYMBOL[[1]]
    weights<-weights[order(abs(weights))]
    gene_ordered <- factor(names(weights), levels = names(weights))
    plot <- data.frame(weight = weights, gene = gene_ordered) %>%
      mutate(Direction = ifelse(weight > 0, "Positive", "Negative")) %>%
      ggplot(aes(x = gene, y = abs(weight), fill = Direction)) + 
      geom_bar(stat = "identity") +
      theme_bw() +
      ylab("Weight") +
      xlab("Genes") +
      ggtitle(paste0(GO_term, "\n"," ", gsub("_", " ", tissue_name))) +  # "\n" solo para selected
      coord_cartesian(ylim = c(0, y_axis_limit)) +  
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 6),
            axis.title.y = element_text(size = 12),  
            legend.text = element_text(size = 11),  
            legend.title = element_text(size = 12)) 
    
    plot_list[[j]] <- plot
  }
  
  combined_plot <- wrap_plots(plot_list, guides = "collect") +
    plot_layout(guides = 'collect') &
    theme(legend.position = "right")
  
  combined_plot_with_label <- ggdraw() +
    draw_plot(combined_plot) +  # The combined boxplot
    draw_label("B", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 15)  # Single "A" label in top-left
  
  # extraemos el summarisedexperiment de uno de los tejidos de interes (que no sea tejido de organos sexuales)
  complete_odds_df<-data.frame()
  
  for (k in 1:nrow(all_pvalues_2))
  {
    tissue_name<-all_pvalues_2$tissue[[k]]
    GO_term<-all_pvalues_2$GO[[k]]
    
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
    conf_int <- confint(model,level = 0.95)  
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
  # Create boxplots for significant GO terms
  boxplot_list <- list()
  color_palette <- c("control" = "#1f77b4",  
                     "AD" = "#ff7f0e")    
  for (l in 1:nrow(all_pvalues_2)) {
    tissue_name <- all_pvalues_2$tissue[[l]]
    GO_term <- all_pvalues_2$GO[[l]]
    
    tissue_score <- score_tissue_list[[tissue_name]]
    
    plot <- data.frame(Expression = as.vector(assay(tissue_score)[GO_term, ]),
                       cases = tissue_score$DiseaseStatus) %>%
      ggplot(aes(x = cases, y = Expression, col = cases)) +
      scale_color_manual(values = color_palette) +  # Apply custom colors
      geom_boxplot() +
      theme_bw() + 
      coord_cartesian(ylim = c(-3, 2)) +
      ylab("GSAS") +
      xlab(NULL) +
      labs(color = "Cases") +
      theme(axis.text.x = element_blank(),
            axis.title.y = element_text(size = 16), 
            legend.position = "right") +
      ggtitle(paste0(GO_term, "\n", gsub("_", " ", tissue_name))) 
    
    boxplot_list[[l]] <- plot
  }
  
  # Combine the boxplots into a grid using patchwork
  combined_boxplot <- wrap_plots(boxplot_list, guides = "collect") +
    plot_layout(guides = 'collect') &
    theme(legend.position = "right", 
          legend.text = element_text(size = 14),  
          legend.title = element_text(size = 16)) 
  
  # Use cowplot to add the label "A" outside the entire plot grid
  combined_boxplot_with_label <- ggdraw() +
    draw_plot(combined_boxplot) +  # The combined boxplot
    draw_label("A", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 15)  # Single "A" label in top-left
  
  # Print the final plot with the label
  print(combined_plot_with_label)
  return(list(
    linear_reg = all_pvalues_2,
    log_reg = complete_odds_df,
    combined_plot_with_label,
    combined_boxplot_with_label))
}



load(paste0(path,  "results/netactivity/se_list_curated_TODO.RData")) # netactivity ejecutar todo
all_tissues_p_value<-pvalue_boxplot(score_tissue_list)



lin_analysis<-all_tissues_p_value[[1]]
reg_analysis<-all_tissues_p_value[[2]]
plot_all<-all_tissues_p_value[[3]]
box_all<-all_tissues_p_value[[4]]

# guardamos la regresion lineal para analizar los genes de estos GOs más adelante 
save(lin_analysis,file=paste0(path,  "results/netactivity/linear_reg.RData")) 

ggsave(paste0(path,"figures/netactivity_images/boxplot_all_tissues.png"), plot = box_all, width = 9, height = 6, dpi = 400)

ggsave(paste0(path, "figures/netactivity_images/plotbar_all_tissues.png"), plot = plot_all, width = 9, height = 6, dpi = 300)


# Ahora representamos los GOs en regresion lineal

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

graph1 <- ggplot(reg_analysis, aes(x = Odds_Ratio, y = reorder(paste(GO, tissue, sep = "_"), Odds_Ratio), fill = GO)) +
  geom_bar(stat = "identity", position = "dodge", aes(x = Odds_Ratio - 1)) +
  
  geom_errorbar(aes(xmin = CI_Lower-1, xmax = CI_Upper-1), width = 0.2, position = position_dodge(0.9)) +
  
  labs(title = "Odds Ratios for GSAS of GO Terms",
       x = "Odds Ratio",
       y = NULL,
       fill = "GO Term") +
  
  scale_fill_manual(values = c("GO:0051481" = "#F8766D",
                               "GO:0002467" = "#619CFF",
                               "GO:0070572" = "#C77CFF", 
                               "GO:0001516" = "#7CAE00")) +
  
  scale_x_continuous(limits = c(-0.5,1), breaks = seq(-0.5, 1, by = 0.25), labels = seq(0.5, 2, by = 0.25)) +
  
  scale_y_discrete(
    limits = rev(unique(paste(reg_analysis$GO, reg_analysis$tissue, sep = "_"))),
    labels = NULL) +
  theme(
    axis.title.x = element_text(size = 10), 
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 12),
  )


hjust <- ifelse(reg_analysis$Odds_Ratio < 1, 1.25, -0.25) 

graph1 + 
  geom_text(
    aes(x = Odds_Ratio - 1, label = signif),
    hjust = hjust, 
    size = 4, 
    na.rm = TRUE, 
    position = position_dodge(1)
  )

ggsave(paste0(path, "figures/netactivity_images/OR_GSAS.png"), plot = graph1, width = 8, height = 3, dpi = 300)

