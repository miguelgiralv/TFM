NetActivity
================
Miguel G.A.
2024-06-24

    ## Loading required package: MatrixGenerics

    ## Warning: package 'MatrixGenerics' was built under R version 4.3.1

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 4.3.3

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 4.3.1

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 4.3.1

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 4.3.1

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 4.3.1

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    ## Warning: package 'tidyverse' was built under R version 4.3.3

    ## Warning: package 'ggplot2' was built under R version 4.3.3

    ## Warning: package 'tidyr' was built under R version 4.3.3

    ## Warning: package 'readr' was built under R version 4.3.3

    ## Warning: package 'purrr' was built under R version 4.3.3

    ## Warning: package 'dplyr' was built under R version 4.3.3

    ## Warning: package 'stringr' was built under R version 4.3.3

    ## Warning: package 'forcats' was built under R version 4.3.1

    ## Warning: package 'lubridate' was built under R version 4.3.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()        masks matrixStats::count()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
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
```

    ## [1] "Adipose_Subcutaneous procesado"
    ## [1] "Adipose_Visceral_Omentum procesado"
    ## [1] "Adrenal_Gland procesado"
    ## [1] "Artery_Aorta procesado"
    ## [1] "Artery_Coronary procesado"
    ## [1] "Artery_Tibial procesado"
    ## [1] "Brain_Amygdala procesado"
    ## [1] "Brain_Anterior_cingulate_cortex_BA24 procesado"
    ## [1] "Brain_Caudate_basal_ganglia procesado"
    ## [1] "Brain_Cerebellar_Hemisphere procesado"
    ## [1] "Brain_Cerebellum procesado"
    ## [1] "Brain_Cortex procesado"
    ## [1] "Brain_Frontal_Cortex_BA9 procesado"
    ## [1] "Brain_Hippocampus procesado"
    ## [1] "Brain_Hypothalamus procesado"
    ## [1] "Brain_Nucleus_accumbens_basal_ganglia procesado"
    ## [1] "Brain_Putamen_basal_ganglia procesado"
    ## [1] "Brain_Spinal_cord_cervical_c-1 procesado"
    ## [1] "Brain_Substantia_nigra procesado"
    ## [1] "Breast_Mammary_Tissue procesado"
    ## [1] "Cells_Cultured_fibroblasts procesado"
    ## [1] "Cells_EBV-transformed_lymphocytes procesado"
    ## [1] "Colon_Sigmoid procesado"
    ## [1] "Colon_Transverse procesado"
    ## [1] "Esophagus_Gastroesophageal_Junction procesado"
    ## [1] "Esophagus_Mucosa procesado"
    ## [1] "Esophagus_Muscularis procesado"
    ## [1] "Heart_Atrial_Appendage procesado"
    ## [1] "Heart_Left_Ventricle procesado"
    ## [1] "Kidney_Cortex procesado"
    ## [1] "Liver procesado"
    ## [1] "Lung procesado"
    ## [1] "Minor_Salivary_Gland procesado"
    ## [1] "Muscle_Skeletal procesado"
    ## [1] "Nerve_Tibial procesado"
    ## [1] "Ovary procesado"
    ## [1] "Pancreas procesado"
    ## [1] "Pituitary procesado"
    ## [1] "Prostate procesado"
    ## [1] "Skin_Not_Sun_Exposed_Suprapubic procesado"
    ## [1] "Skin_Sun_Exposed_Lower_leg procesado"
    ## [1] "Small_Intestine_Terminal_Ileum procesado"
    ## [1] "Spleen procesado"
    ## [1] "Stomach procesado"
    ## [1] "Testis procesado"
    ## [1] "Thyroid procesado"
    ## [1] "Uterus procesado"
    ## [1] "Vagina procesado"
    ## [1] "Whole_Blood procesado"

``` r
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
```

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2424 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Adipose_Subcutaneous computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2466 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Adipose_Visceral_Omentum computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2926 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Adrenal_Gland computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2492 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Artery_Aorta computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2848 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Artery_Coronary computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2420 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Artery_Tibial computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3396 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Amygdala computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3108 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Anterior_cingulate_cortex_BA24 computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2883 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Caudate_basal_ganglia computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2945 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Cerebellar_Hemisphere computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2862 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Cerebellum computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2806 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Cortex computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2852 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Frontal_Cortex_BA9 computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3107 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Hippocampus computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3060 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Hypothalamus computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2915 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Nucleus_accumbens_basal_ganglia computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2994 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Putamen_basal_ganglia computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3224 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Spinal_cord_cervical_c-1 computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3441 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Brain_Substantia_nigra computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2566 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Breast_Mammary_Tissue computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2568 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Cells_Cultured_fibroblasts computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3277 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Cells_EBV-transformed_lymphocytes computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2644 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Colon_Sigmoid computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2580 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Colon_Transverse computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2593 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Esophagus_Gastroesophageal_Junction computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2421 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Esophagus_Mucosa computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2415 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Esophagus_Muscularis computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2686 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Heart_Atrial_Appendage computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2915 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Heart_Left_Ventricle computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 4076 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Kidney_Cortex computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3205 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Liver computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2390 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Lung computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2904 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Minor_Salivary_Gland computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2725 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Muscle_Skeletal computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2227 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Nerve_Tibial computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2911 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Ovary computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2754 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Pancreas computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2714 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Pituitary computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2773 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Prostate computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2326 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Skin_Not_Sun_Exposed_Suprapubic computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2253 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Skin_Sun_Exposed_Lower_leg computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2852 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Small_Intestine_Terminal_Ileum computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2853 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Spleen computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2699 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Stomach computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2067 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Testis computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 2265 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Thyroid computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3183 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Uterus computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3335 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Vagina computado"

    ## Warning in prepareSummarizedExperiment(se_tissue, "gtex_gokegg"): 3031 genes
    ## present in the model not found in input data. The expression of all samples
    ## will be set to 0.

    ## [1] "Whole_Blood computado"

``` r
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
```

    ## Adipose_Subcutaneous filas con valor de cero: 0

    ## Adipose_Subcutaneous procesado, filas restantes: 1518

    ## Adipose_Visceral_Omentum filas con valor de cero: 0

    ## Adipose_Visceral_Omentum procesado, filas restantes: 1518

    ## Adrenal_Gland filas con valor de cero: 0

    ## Adrenal_Gland procesado, filas restantes: 1518

    ## Artery_Aorta filas con valor de cero: 0

    ## Artery_Aorta procesado, filas restantes: 1518

    ## Artery_Coronary filas con valor de cero: 0

    ## Artery_Coronary procesado, filas restantes: 1518

    ## Artery_Tibial filas con valor de cero: 0

    ## Artery_Tibial procesado, filas restantes: 1518

    ## Brain_Amygdala filas con valor de cero: 0

    ## Brain_Amygdala procesado, filas restantes: 1518

    ## Brain_Anterior_cingulate_cortex_BA24 filas con valor de cero: 0

    ## Brain_Anterior_cingulate_cortex_BA24 procesado, filas restantes: 1518

    ## Brain_Caudate_basal_ganglia filas con valor de cero: 0

    ## Brain_Caudate_basal_ganglia procesado, filas restantes: 1518

    ## Brain_Cerebellar_Hemisphere filas con valor de cero: 0

    ## Brain_Cerebellar_Hemisphere procesado, filas restantes: 1518

    ## Brain_Cerebellum filas con valor de cero: 0

    ## Brain_Cerebellum procesado, filas restantes: 1518

    ## Brain_Cortex filas con valor de cero: 0

    ## Brain_Cortex procesado, filas restantes: 1518

    ## Brain_Frontal_Cortex_BA9 filas con valor de cero: 0

    ## Brain_Frontal_Cortex_BA9 procesado, filas restantes: 1518

    ## Brain_Hippocampus filas con valor de cero: 0

    ## Brain_Hippocampus procesado, filas restantes: 1518

    ## Brain_Hypothalamus filas con valor de cero: 0

    ## Brain_Hypothalamus procesado, filas restantes: 1518

    ## Brain_Nucleus_accumbens_basal_ganglia filas con valor de cero: 1

    ## Brain_Nucleus_accumbens_basal_ganglia procesado, filas restantes: 1517

    ## Brain_Putamen_basal_ganglia filas con valor de cero: 0

    ## Brain_Putamen_basal_ganglia procesado, filas restantes: 1518

    ## Brain_Spinal_cord_cervical_c-1 filas con valor de cero: 0

    ## Brain_Spinal_cord_cervical_c-1 procesado, filas restantes: 1518

    ## Brain_Substantia_nigra filas con valor de cero: 1

    ## Brain_Substantia_nigra procesado, filas restantes: 1517

    ## Breast_Mammary_Tissue filas con valor de cero: 0

    ## Breast_Mammary_Tissue procesado, filas restantes: 1518

    ## Cells_Cultured_fibroblasts filas con valor de cero: 0

    ## Cells_Cultured_fibroblasts procesado, filas restantes: 1518

    ## Cells_EBV-transformed_lymphocytes filas con valor de cero: 0

    ## Cells_EBV-transformed_lymphocytes procesado, filas restantes: 1518

    ## Colon_Sigmoid filas con valor de cero: 0

    ## Colon_Sigmoid procesado, filas restantes: 1518

    ## Colon_Transverse filas con valor de cero: 0

    ## Colon_Transverse procesado, filas restantes: 1518

    ## Esophagus_Gastroesophageal_Junction filas con valor de cero: 0

    ## Esophagus_Gastroesophageal_Junction procesado, filas restantes: 1518

    ## Esophagus_Mucosa filas con valor de cero: 0

    ## Esophagus_Mucosa procesado, filas restantes: 1518

    ## Esophagus_Muscularis filas con valor de cero: 0

    ## Esophagus_Muscularis procesado, filas restantes: 1518

    ## Heart_Atrial_Appendage filas con valor de cero: 0

    ## Heart_Atrial_Appendage procesado, filas restantes: 1518

    ## Heart_Left_Ventricle filas con valor de cero: 0

    ## Heart_Left_Ventricle procesado, filas restantes: 1518

    ## Kidney_Cortex filas con valor de cero: 1

    ## Kidney_Cortex procesado, filas restantes: 1517

    ## Liver filas con valor de cero: 1

    ## Liver procesado, filas restantes: 1517

    ## Lung filas con valor de cero: 1

    ## Lung procesado, filas restantes: 1517

    ## Minor_Salivary_Gland filas con valor de cero: 0

    ## Minor_Salivary_Gland procesado, filas restantes: 1518

    ## Muscle_Skeletal filas con valor de cero: 1

    ## Muscle_Skeletal procesado, filas restantes: 1517

    ## Nerve_Tibial filas con valor de cero: 0

    ## Nerve_Tibial procesado, filas restantes: 1518

    ## Ovary filas con valor de cero: 0

    ## Ovary procesado, filas restantes: 1518

    ## Pancreas filas con valor de cero: 0

    ## Pancreas procesado, filas restantes: 1518

    ## Pituitary filas con valor de cero: 0

    ## Pituitary procesado, filas restantes: 1518

    ## Prostate filas con valor de cero: 0

    ## Prostate procesado, filas restantes: 1518

    ## Skin_Not_Sun_Exposed_Suprapubic filas con valor de cero: 0

    ## Skin_Not_Sun_Exposed_Suprapubic procesado, filas restantes: 1518

    ## Skin_Sun_Exposed_Lower_leg filas con valor de cero: 0

    ## Skin_Sun_Exposed_Lower_leg procesado, filas restantes: 1518

    ## Small_Intestine_Terminal_Ileum filas con valor de cero: 0

    ## Small_Intestine_Terminal_Ileum procesado, filas restantes: 1518

    ## Spleen filas con valor de cero: 0

    ## Spleen procesado, filas restantes: 1518

    ## Stomach filas con valor de cero: 0

    ## Stomach procesado, filas restantes: 1518

    ## Testis filas con valor de cero: 0

    ## Testis procesado, filas restantes: 1518

    ## Thyroid filas con valor de cero: 0

    ## Thyroid procesado, filas restantes: 1518

    ## Uterus filas con valor de cero: 0

    ## Uterus procesado, filas restantes: 1518

    ## Vagina filas con valor de cero: 2

    ## Vagina procesado, filas restantes: 1516

    ## Whole_Blood filas con valor de cero: 1

    ## Whole_Blood procesado, filas restantes: 1517

``` r
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
```

    ## [1] "Adipose_Subcutaneous computado"
    ## [1] "Adipose_Visceral_Omentum computado"
    ## [1] "Adrenal_Gland computado"
    ## [1] "Artery_Aorta computado"
    ## [1] "Artery_Coronary computado"
    ## [1] "Artery_Tibial computado"
    ## [1] "Brain_Amygdala computado"
    ## [1] "Brain_Anterior_cingulate_cortex_BA24 computado"
    ## [1] "Brain_Caudate_basal_ganglia computado"
    ## [1] "Brain_Cerebellar_Hemisphere computado"
    ## [1] "Brain_Cerebellum computado"
    ## [1] "Brain_Cortex computado"
    ## [1] "Brain_Frontal_Cortex_BA9 computado"
    ## [1] "Brain_Hippocampus computado"
    ## [1] "Brain_Hypothalamus computado"
    ## [1] "Brain_Nucleus_accumbens_basal_ganglia computado"
    ## [1] "Brain_Putamen_basal_ganglia computado"
    ## [1] "Brain_Spinal_cord_cervical_c-1 computado"
    ## [1] "Brain_Substantia_nigra computado"
    ## [1] "Breast_Mammary_Tissue computado"
    ## [1] "Cells_Cultured_fibroblasts computado"
    ## [1] "Cells_EBV-transformed_lymphocytes computado"
    ## [1] "Colon_Sigmoid computado"
    ## [1] "Colon_Transverse computado"
    ## [1] "Esophagus_Gastroesophageal_Junction computado"
    ## [1] "Esophagus_Mucosa computado"
    ## [1] "Esophagus_Muscularis computado"
    ## [1] "Heart_Atrial_Appendage computado"
    ## [1] "Heart_Left_Ventricle computado"
    ## [1] "Kidney_Cortex computado"
    ## [1] "Liver computado"
    ## [1] "Lung computado"
    ## [1] "Minor_Salivary_Gland computado"
    ## [1] "Muscle_Skeletal computado"
    ## [1] "Nerve_Tibial computado"
    ## [1] "Ovary computado"
    ## [1] "Pancreas computado"
    ## [1] "Pituitary computado"
    ## [1] "Prostate computado"
    ## [1] "Skin_Not_Sun_Exposed_Suprapubic computado"
    ## [1] "Skin_Sun_Exposed_Lower_leg computado"
    ## [1] "Small_Intestine_Terminal_Ileum computado"
    ## [1] "Spleen computado"
    ## [1] "Stomach computado"
    ## [1] "Testis computado"
    ## [1] "Thyroid computado"
    ## [1] "Uterus computado"
    ## [1] "Vagina computado"
    ## [1] "Whole_Blood computado"

``` r
all_pvalues$Bonferroni <- p.adjust(all_pvalues$DiseaseStatusAD, method = "bonferroni")  
all_pvalues$FDR <- p.adjust(all_pvalues$DiseaseStatusAD, method = "fdr")  

all_pvalues_Bonferroni<-all_pvalues %>% select("DiseaseStatusAD", "tissue","Bonferroni")
all_pvalues_Bonferroni<-all_pvalues_Bonferroni[order(all_pvalues_Bonferroni$Bonferroni), ]

all_pvalues_FDR<-all_pvalues %>% select("DiseaseStatusAD", "tissue","FDR")
all_pvalues_FDR <- all_pvalues_FDR[order(all_pvalues_FDR$FDR), ]

head(all_pvalues_Bonferroni,20)
```

    ##              DiseaseStatusAD               tissue Bonferroni
    ## GO:000246726    6.388980e-06 Esophagus_Muscularis  0.4751676
    ## GO:0006734      6.141044e-02 Adipose_Subcutaneous  1.0000000
    ## GO:0008340      5.350899e-02 Adipose_Subcutaneous  1.0000000
    ## GO:0014854      5.171780e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0019081      2.703077e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0019835      1.555432e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0031987      2.512516e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0042753      2.362277e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0042754      4.320771e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0043476      2.940238e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0044793      3.106875e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0044794      3.594915e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0045056      4.306501e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0045780      5.311531e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0046755      2.265199e-02 Adipose_Subcutaneous  1.0000000
    ## GO:0046794      2.407591e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0048532      8.892604e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0051775      3.583143e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0060033      6.597771e-01 Adipose_Subcutaneous  1.0000000
    ## GO:0060352      4.565938e-01 Adipose_Subcutaneous  1.0000000

``` r
head(all_pvalues_FDR,20)
```

    ##              DiseaseStatusAD                              tissue       FDR
    ## GO:000246726    6.388980e-06                Esophagus_Muscularis 0.4751676
    ## GO:19000342     3.841028e-05                       Adrenal_Gland 0.5739564
    ## GO:005148113    2.476395e-05                   Brain_Hippocampus 0.5739564
    ## GO:001975525    3.858634e-05                    Esophagus_Mucosa 0.5739564
    ## GO:007057228    3.584486e-05                Heart_Left_Ventricle 0.5739564
    ## GO:0072576      1.768841e-04                Adipose_Subcutaneous 0.5965195
    ## GO:00705721     1.160306e-04            Adipose_Visceral_Omentum 0.5965195
    ## GO:00218713     6.803708e-05                        Artery_Aorta 0.5965195
    ## GO:00725763     1.414772e-04                        Artery_Aorta 0.5965195
    ## GO:00424034     1.424671e-04                     Artery_Coronary 0.5965195
    ## GO:00725764     9.571409e-05                     Artery_Coronary 0.5965195
    ## GO:00601924     1.173823e-04                     Artery_Coronary 0.5965195
    ## GO:004240311    1.627330e-04                        Brain_Cortex 0.5965195
    ## GO:007154920    1.100644e-04          Cells_Cultured_fibroblasts 0.5965195
    ## GO:001580124    8.250311e-05 Esophagus_Gastroesophageal_Junction 0.5965195
    ## GO:007257626    1.040098e-04                Esophagus_Muscularis 0.5965195
    ## GO:001975527    1.528368e-04              Heart_Atrial_Appendage 0.5965195
    ## GO:000151628    5.736018e-05                Heart_Left_Ventricle 0.5965195
    ## GO:000246729    1.622673e-04                       Kidney_Cortex 0.5965195
    ## GO:003332732    1.568137e-04                Minor_Salivary_Gland 0.5965195

``` r
write.csv(all_pvalues_Bonferroni, file = paste0(path,"results/netactivity/", "bonferroni.txt"), row.names = TRUE)
write.csv(all_pvalues_FDR, file = paste0(path,"results/netactivity/","FDR.txt"), row.names = TRUE)
```
