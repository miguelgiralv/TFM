NetActivity
================
Miguel G.A.
2024-06-24

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

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

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

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

    ## Warning: package 'SparseArray' was built under R version 4.3.1

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.3.1

    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand
    ## 
    ## Loading required package: S4Arrays

    ## Warning: package 'S4Arrays' was built under R version 4.3.1

    ## Loading required package: abind
    ## 
    ## Attaching package: 'S4Arrays'
    ## 
    ## The following object is masked from 'package:abind':
    ## 
    ##     abind
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     rowsum
    ## 
    ## Registered S3 method overwritten by 'SparseArray':
    ##   method           from        
    ##   rowsum.dgCMatrix DelayedArray

Cargamos la expresion en el tejido adiposo subcutaneo como un df en R:

``` r
df_adipose_en <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_en/processed/processed_en_Adipose_Subcutaneous_predict.txt", header = TRUE, sep = "\t")


dim(df_adipose_en)
```

    ## [1]  878 8635

Procesamos las filas para que en los individuos no aparezca el gse y
tambien cambiamos los nombres de los genes para que no aparezca la
version de ensembl

``` r
individuals <- df_adipose_en[, 1]
individuals <- sub("^GSE33528_", "", individuals)

df_adipose_en[, 1] <- individuals

colnames(df_adipose_en) <- sub("\\.\\d+$", "", colnames(df_adipose_en))
individuals_cases <- read.csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", header = TRUE, sep = "\t")
```

Fusionamos las dos matrices de individuos con metadatos y creamos una
nueva que contenga solo los datos del experimento y su enfermedad

``` r
merged_individuals <- merge(
  x = data.frame(Individual.ID = individuals), 
  y = individuals_cases[, c("Individual.ID", "Disease.status")],
  by.x = "Individual.ID",
  by.y = "Individual.ID",
  all.x = TRUE  
)
```

Ahora transponemos nuestra matriz de expresion de forma que para cada
fila haya un gen y para cada individuo haya una columna

``` r
gene_columns <- df_adipose_en[, !names(df_adipose_en) %in% "IID", drop = FALSE]

transposed_genes <- t(as.matrix(gene_columns))

df_adipose_en_trans <- as.data.frame(transposed_genes)


colnames(df_adipose_en_trans)<-individuals


dim(df_adipose_en_trans)
```

    ## [1] 8634  878

Finalmente añadimos row_data y col_data y cargamos todos los datos para
cargar el objeto SummarizedExperiment

``` r
genes <- df_adipose_en[, !names(df_adipose_en) %in% "IID"]


row_data <- DataFrame(
  Gene = rownames(df_adipose_en_trans)
)

col_data <- data.frame(
  Sample = merged_individuals$Individual.ID,
  DiseaseStatus = merged_individuals$Disease.status
)



se_adipose_en <- SummarizedExperiment(
  assays = list(expression = as.matrix(df_adipose_en_trans)), # lo guardamos como una matriz y no un df
  colData = col_data,
  rowData = row_data
)

str(se_adipose_en)
```

    ## Formal class 'SummarizedExperiment' [package "SummarizedExperiment"] with 5 slots
    ##   ..@ colData        :Formal class 'DFrame' [package "S4Vectors"] with 6 slots
    ##   .. .. ..@ rownames       : chr [1:878] "GSM894503" "GSM894504" "GSM894505" "GSM894506" ...
    ##   .. .. ..@ nrows          : int 878
    ##   .. .. ..@ elementType    : chr "ANY"
    ##   .. .. ..@ elementMetadata: NULL
    ##   .. .. ..@ metadata       : list()
    ##   .. .. ..@ listData       :List of 2
    ##   .. .. .. ..$ Sample       : chr [1:878] "GSM894503" "GSM894504" "GSM894505" "GSM894506" ...
    ##   .. .. .. ..$ DiseaseStatus: chr [1:878] "control" "control" "AD" "control" ...
    ##   ..@ assays         :Formal class 'SimpleAssays' [package "SummarizedExperiment"] with 1 slot
    ##   .. .. ..@ data:Formal class 'SimpleList' [package "S4Vectors"] with 4 slots
    ##   .. .. .. .. ..@ listData       :List of 1
    ##   .. .. .. .. .. ..$ expression: num [1:8634, 1:878] 0.0549 -0.2539 -0.0368 0.0765 0.1023 ...
    ##   .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. .. .. .. .. .. ..$ : chr [1:8634] "ENSG00000237491" "ENSG00000225880" "ENSG00000248333" "ENSG00000272512" ...
    ##   .. .. .. .. .. .. .. ..$ : chr [1:878] "GSM894503" "GSM894504" "GSM894505" "GSM894506" ...
    ##   .. .. .. .. ..@ elementType    : chr "ANY"
    ##   .. .. .. .. ..@ elementMetadata: NULL
    ##   .. .. .. .. ..@ metadata       : list()
    ##   ..@ NAMES          : chr [1:8634] "ENSG00000237491" "ENSG00000225880" "ENSG00000248333" "ENSG00000272512" ...
    ##   ..@ elementMetadata:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
    ##   .. .. ..@ rownames       : NULL
    ##   .. .. ..@ nrows          : int 8634
    ##   .. .. ..@ elementType    : chr "ANY"
    ##   .. .. ..@ elementMetadata: NULL
    ##   .. .. ..@ metadata       : list()
    ##   .. .. ..@ listData       :List of 1
    ##   .. .. .. ..$ Gene: chr [1:8634] "ENSG00000237491" "ENSG00000225880" "ENSG00000248333" "ENSG00000272512" ...
    ##   ..@ metadata       : list()

Ahora podemos ya aplicar netactivity

``` r
out_adipose_en <- prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg")
```

    ## Warning in prepareSummarizedExperiment(se_adipose_en, "gtex_gokegg"): 5031
    ## genes present in the model not found in input data. The expression of all
    ## samples will be set to 0.

``` r
scores_adipose_en <- computeGeneSetScores(out_adipose_en, "gtex_gokegg")
```

``` r
mod_adipose_en <- model.matrix(~ DiseaseStatus , colData(scores_adipose_en))
fit <- lmFit(assay(scores_adipose_en), mod_adipose_en) %>% eBayes()
```

    ## Warning: Zero sample variances detected, have been offset away from zero

``` r
topTab <- topTable(fit, coef = 1:2, n = Inf)
head(topTab)
```

    ##            X.Intercept. DiseaseStatuscontrol       AveExpr        F     P.Value
    ## GO:0090200   0.04412589          -0.09202501  4.168868e-19 5.852282 0.002986537
    ## GO:0046322  -0.02980884           0.06216665 -1.678019e-17 4.751058 0.008865852
    ## GO:0071800  -0.02342122           0.04884520 -1.728302e-17 4.641852 0.009877483
    ## GO:0035815  -0.03049967           0.06360740  4.508700e-18 4.608886 0.010205022
    ## GO:0008334   0.03364799          -0.07017323 -1.186942e-17 4.302707 0.013818170
    ## GO:0071549   0.02451813          -0.05113282 -8.811921e-19 4.267400 0.014309878
    ##            adj.P.Val
    ## GO:0090200         1
    ## GO:0046322         1
    ## GO:0071800         1
    ## GO:0035815         1
    ## GO:0008334         1
    ## GO:0071549         1
