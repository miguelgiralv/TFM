library(tidyverse)
library(tibble)
library(dplyr)
library(caret)
library(grid)
library(pROC)

# el path al repositorio

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" 

# el nuevo creado
load(paste0(path, "results/netactivity/Metascore_TODO.RData"))
load(paste0(path, "results/netactivity/Metascore_selected.RData"))
load(paste0(path,  "results/PRS.RData"))
load(paste0(path,  "results/PRS.RData"))

Metascore_test_list_all[["fold_1"]][-3]
Metascore_test_list_selected[["fold_1"]][-3]

## PRIMERO REPRESENTAMOS Y COMPARAMOS EL RENDIMIENTO DE todos los MODELOs DE FORMA INDIVIDUAL
# PARA EL TODO:
metrics_summary_all<-Metascore_test_list_all[["Performance"]]
metrics_summary_selected<-Metascore_test_list_selected[["Performance"]]

metrics_summary_all$Dataset <- "All tissues"
metrics_summary_selected$Dataset <- "Selected tissues"

# guardamos las tablas:
metrics_summary_all_long <- metrics_summary_all %>%
  filter(Fold != "Average") %>%
  pivot_longer(cols = -c(Fold, Dataset), names_to = "Metric", values_to = "Value")

metrics_summary_selected_long <- metrics_summary_selected %>%
  filter(Fold != "Average") %>%
  pivot_longer(cols = -c(Fold, Dataset), names_to = "Metric", values_to = "Value")

# sacamos las metricas de los PRS
prs_results<-prs_results[-5,] # quitamos el f score

PRS_results_general_OR<-PRS_results_general[[1]]
PRS_results_general_metrics<-PRS_results_general[[2]]
PRS_results_general_metrics<-PRS_results_general_metrics[-5,]

PRS_results_hispana_OR<-PRS_results_hispana[[1]]
PRS_results_hispana_metrics<-PRS_results_hispana[[2]]
PRS_results_hispana_metrics<-PRS_results_hispana_metrics[-5,]

# Cogemos las metricas de los PRS:
PRS_results_general_metrics <- PRS_results_general_metrics %>%
  mutate(Dataset = "PGS004146")

PRS_results_hispana_metrics <- PRS_results_hispana_metrics %>%
  mutate(Dataset = "PGS000054")

metrics_summary_all_long <- metrics_summary_all %>%
  filter(Fold != "Average") %>%
  pivot_longer(cols = -c(Fold, Dataset), names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "All tissues")

metrics_summary_selected_long <- metrics_summary_selected %>%
  filter(Fold != "Average") %>%
  pivot_longer(cols = -c(Fold, Dataset), names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "Selected tissues")

# Combinamos 
metrics_combined2 <- bind_rows(
  metrics_summary_all_long,
  metrics_summary_selected_long,
  PRS_results_general_metrics,
  PRS_results_hispana_metrics
)

# BOXPLOT REPRESENTACION:
boxplot_compared<-ggplot(metrics_combined2, aes(x = Metric, y = Value, fill = Dataset)) +
  # metemos en el boxplot todo excepto los prs
  geom_boxplot(data = metrics_combined2 %>% filter(Dataset != "PGS004146" & Dataset != "PGS000054"), 
               position = position_dodge(width = 0.8)) + 
  
  # lineas discontinuas para los prs
  geom_segment(data = metrics_combined2 %>% filter(Dataset %in% c("PGS004146", "PGS000054")), 
               aes(x = as.numeric(as.factor(Metric)) - 0.43, 
                   xend = as.numeric(as.factor(Metric)) + 0.43, 
                   y = Value, yend = Value, color = Dataset), 
               linetype = "dotted", size = 1.5) + 
  labs(
    title = "Comparison of Performance Metrics Across Folds",
    x = "Metrics",
    y = "Value",
    fill = "Models",  
    color = "PRS"  
  ) +
  coord_cartesian(ylim = c(0.15, 0.85)) + 
  scale_fill_manual(
    values = c("All tissues" = "#F8766D",
               "Selected tissues" = "#00BFC4"),
    breaks = c("All tissues",
               "Selected tissues")
  ) +
  scale_color_manual(
    values = c("PGS004146" = "#C77CFF", 
               "PGS000054" = "#7CAE00"),
    
    breaks = c("PGS004146", "PGS000054")) +
  theme(
  axis.title.y = element_text(size = 14),  
  axis.title.x = element_text(size = 14),  #
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  plot.title = element_text(size = 16),
)

ggsave(paste0(path, "figures/metascore/boxplot_PRS_metascores.png"), plot = boxplot_compared, width = 9, height = 6, dpi = 400)

# ahora representamos los folds y los odds ratio

comparison_list<-list()

for (fold in 1:5)
{
  Metascore_test_list_all_OR<-Metascore_test_list_all[[fold]][["odds_ratio"]]
  rownames(Metascore_test_list_all_OR)<-c("Intercept_all","Metascore_all")
  Metascore_test_list_selected_OR<-Metascore_test_list_selected[[fold]][["odds_ratio"]]
  rownames(Metascore_test_list_selected_OR)<-c("Intercept_selected","Metascore_selected")
  combined_odds_ratio<-rbind(Metascore_test_list_all_OR,Metascore_test_list_selected_OR)
  comparison_list[[paste0("fold_", fold)]]<-combined_odds_ratio
}

# meter los asteriscos tb  *<0.1 **<0.05 ***<0.01
main_data_total_df<-data.frame()
barplot_list<-list()

for (fold in names(comparison_list))
{
  estimates_df<-comparison_list[[fold]]
  estimates_df<-as.data.frame(estimates_df)
  colnames(estimates_df)[7]<-"p_value"
  estimates_df$signif <- ifelse(
    estimates_df$p_value > 0.1, NA, 
    ifelse(
      estimates_df$p_value < 0.1 & estimates_df$p_value > 0.05, "*", 
      ifelse(
        estimates_df$p_value < 0.05 & estimates_df$p_value > 0.01, "**", 
        ifelse(
          estimates_df$p_value < 0.01, "***", NA
        )
      )
    )
  )
  estimates_df$fold<-fold
  estimates_df$model<-rownames(estimates_df)
  estimates_df$complete_name<-paste0(rownames(estimates_df),"_", fold)
  rownames(estimates_df)<-NULL
  main_data_total_df<-rbind(main_data_total_df,estimates_df)
}


all_OR<-subset(filtered_data, model == "Metascore_all")
all_sel<-subset(filtered_data, model == "Metascore_selected")

mean(all_sel$Odds_Ratio)
mean(all_OR$Odds_Ratio)


filtered_data$model <- factor(filtered_data$model, 
                             levels = c("Metascore_all", "Metascore_selected"))
filtered_data

colnames(PRS_results_general_OR)[7]<-"p_value"

PRS_results_general_OR$signif <- ifelse(
  PRS_results_general_OR$p_value > 0.1, NA, 
  ifelse(
    PRS_results_general_OR$p_value < 0.1 & PRS_results_general_OR$p_value > 0.05, "*", 
    ifelse(
      PRS_results_general_OR$p_value < 0.05 & PRS_results_general_OR$p_value > 0.01, "**", 
      ifelse(
        PRS_results_general_OR$p_value < 0.01, "***", NA
      )
    )
  )
)

PRS_results_general_OR<-PRS_results_general_OR[-1,]
colnames(PRS_results_hispana_OR)[7]<-"p_value"
PRS_results_hispana_OR$signif <- ifelse(
  PRS_results_hispana_OR$p_value > 0.1, NA, 
  ifelse(
    PRS_results_hispana_OR$p_value < 0.1 & PRS_results_hispana_OR$p_value > 0.05, "*", 
    ifelse(
      PRS_results_hispana_OR$p_value < 0.05 & PRS_results_hispana_OR$p_value > 0.01, "**", 
      ifelse(
        PRS_results_hispana_OR$p_value < 0.01, "***", NA
      )
    )
  )
)

PRS_results_hispana_OR<-PRS_results_hispana_OR[-1,]
PRS_results_general_OR <- cbind(PRS_results_general_OR, fold = NA, model = "PGS004146", complete_name = "PGS004146")
PRS_results_hispana_OR <- cbind(PRS_results_hispana_OR, fold = NA, model = "PGS000054", complete_name =  "PGS000054")

combined_data <- rbind(
  filtered_data %>% mutate(Source = "Filtered Data"),
  PRS_results_general_OR %>% mutate(Source = "PRS Results General"),
  PRS_results_hispana_OR %>% mutate(Source = "PRS Results Hispana")
)

combined_data <- combined_data %>% filter(!is.na(model))
graph1 <-   ggplot(combined_data, aes(x = Odds_Ratio, y = complete_name, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", aes(x = Odds_Ratio - 1)) +
  
  # las barras del IC
  geom_errorbar(aes(xmin = CI_Lower - 1, xmax = CI_Upper - 1), width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Odds Ratios of Metascores",
       x = "Odds Ratio",
       y = "Models",
       fill = "Models") +
  scale_fill_manual(values = c("Metascore_all" = "#F8766D",
                               "Metascore_selected" = "#619CFF",
                               "PGS004146" = "#C77CFF", 
                               "PGS000054" = "#7CAE00"),
                    labels = c("Metascore_all" = "All tissues",
                               "Metascore_selected" = "Selected tissues",
                               "PGS004146" = "PGS004146", 
                               "PGS000054" = "PGS000054")) +
  scale_x_continuous(limits = c(-0.5,1), breaks = seq(-0.5, 1, by = 0.25), labels = seq(0.5, 2, by = 0.25)) +
  scale_y_discrete(
    limits = c(
      "Metascore_all_fold_1", "Metascore_all_fold_2", "Metascore_all_fold_3", "Metascore_all_fold_4", "Metascore_all_fold_5",
      "Metascore_selected_fold_1", "Metascore_selected_fold_2", "Metascore_selected_fold_3", "Metascore_selected_fold_4", "Metascore_selected_fold_5",
      "PGS004146", "PGS000054"
    ),
    labels = c(rep(c("Split 1", "Split 2", "Split 3", "Split 4", "Split 5"),2),"PGS004146", "PGS000054")) +
  
  theme(
    axis.title.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16),
  )

hjust <- ifelse(combined_data$Odds_Ratio < 1, 1.25, -0.25) #para los valores menores de 1 le ponemos un hjust hacia la izqda (1.5)
graph1<-graph1 + 
  geom_text(
    aes(x = Odds_Ratio - 1, label = signif),
    hjust = hjust, 
    size = 4, 
    na.rm = TRUE, 
    position = position_dodge(1)
    )
ggsave(paste0(path, "figures/metascore/barplot_PRS_metascores.png"), plot = graph1, width = 9.5, height = 6, dpi = 400)
