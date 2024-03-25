Explore_GSE
================
Miguel G.A.
2024-03-25

## R Markdown

Analizamos el GPL para obtener datos sobre el experimento:

``` r
gpl <- getGEO("GPL14932")

str(gpl)
```

    ## Formal class 'GPL' [package "GEOquery"] with 2 slots
    ##   ..@ dataTable:Formal class 'GEODataTable' [package "GEOquery"] with 2 slots
    ##   .. .. ..@ columns:'data.frame':    16 obs. of  2 variables:
    ##   .. .. .. ..$ Column     : chr [1:16] "ID" "" "" "SNP_ID" ...
    ##   .. .. .. ..$ Description: chr [1:16] "IlmnID" "Name =" "IlmnStrand =" "ncbi SNP accession LINK_PRE:\"http://www.ncbi.nlm.nih.gov/snp/?term=\"" ...
    ##   .. .. ..@ table  :'data.frame':    662841 obs. of  16 variables:
    ##   .. .. .. ..$ ID              : chr [1:662841] "MitoA10045G-gi13273284_B_R_IFB1141652022:" "MitoA10551G-gi13273286_T_F_IFB1141639111:" "MitoA11252G-gi13273288_B_R_IFB1141661584:" "MitoA11468G-gi13273289_B_R_IFB1141719645:" ...
    ##   .. .. .. ..$ Name            : chr [1:662841] "MitoA10045G" "MitoA10551G" "MitoA11252G" "MitoA11468G" ...
    ##   .. .. .. ..$ IlmnStrand      : chr [1:662841] "Bot" "Top" "Bot" "Bot" ...
    ##   .. .. .. ..$ SNP_ID          : chr [1:662841] "[T/C]" "[A/G]" "[T/C]" "[T/C]" ...
    ##   .. .. .. ..$ AddressA_ID     : int [1:662841] 904860368 905860131 903990168 901990681 904780630 906590470 901740201 901470575 905910687 906940170 ...
    ##   .. .. .. ..$ AlleleA_ProbeSeq: logi [1:662841] NA NA NA NA NA NA ...
    ##   .. .. .. ..$ AddressB_ID     : int [1:662841] 0 0 0 0 0 0 0 0 0 0 ...
    ##   .. .. .. ..$ AlleleB_ProbeSeq: logi [1:662841] NA NA NA NA NA NA ...
    ##   .. .. .. ..$ Chr             : chr [1:662841] "M" "M" "M" "M" ...
    ##   .. .. .. ..$ MapInfo         : int [1:662841] 10045 10551 11252 11468 11813 12309 13106 13264 13781 14234 ...
    ##   .. .. .. ..$ Ploidy          : chr [1:662841] "1000" "1000" "1000" "1000" ...
    ##   .. .. .. ..$ Species         : chr [1:662841] "Homo sapiens" "Homo sapiens" "Homo sapiens" "Homo sapiens" ...
    ##   .. .. .. ..$ CustomerStrand  : chr [1:662841] "Top" "Top" "Top" "Top" ...
    ##   .. .. .. ..$ IlmnStrand      : chr [1:662841] "Bot" "Top" "Bot" "Bot" ...
    ##   .. .. .. ..$ IllumicodeSeq   : logi [1:662841] NA NA NA NA NA NA ...
    ##   .. .. .. ..$ TopGenomicSeq   : logi [1:662841] NA NA NA NA NA NA ...
    ##   ..@ header   :List of 25
    ##   .. ..$ contact_address        : chr "6 Queen’s Park Crescent West"
    ##   .. ..$ contact_city           : chr "Toronto"
    ##   .. ..$ contact_country        : chr "Canada"
    ##   .. ..$ contact_email          : chr "ekaterina.rogaeva@utoronto.ca"
    ##   .. ..$ contact_fax            : chr "(416)-978-1878"
    ##   .. ..$ contact_institute      : chr "Tanz Centre for Neurodegenerative Diseases"
    ##   .. ..$ contact_name           : chr "Ekaterina,,Rogaeva"
    ##   .. ..$ contact_phone          : chr "(416) 946-7927"
    ##   .. ..$ contact_zip/postal_code: chr "M5S 3H2"
    ##   .. ..$ data_row_count         : chr "662841"
    ##   .. ..$ description            : chr [1:3] "Descriptor File Name(s),HumanHap650Y_v2-0_GenotypingBC_11237679_A.bpm" "Assay Format,Infinium II" "SNP Count,662841"
    ##   .. ..$ distribution           : chr "custom-commercial"
    ##   .. ..$ geo_accession          : chr "GPL14932"
    ##   .. ..$ last_update_date       : chr "May 01 2012"
    ##   .. ..$ manufacture_protocol   : chr "See manufacturer's website"
    ##   .. ..$ manufacturer           : chr "Illumina, Inc."
    ##   .. ..$ organism               : chr "Homo sapiens"
    ##   .. ..$ sample_id              : chr [1:1215] "GSM894502" "GSM894503" "GSM894504" "GSM894505" ...
    ##   .. ..$ series_id              : chr "GSE33528"
    ##   .. ..$ status                 : chr "Public on May 01 2012"
    ##   .. ..$ submission_date        : chr "Nov 26 2011"
    ##   .. ..$ supplementary_file     : chr "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14932/suppl/GPL14932_Mahdi_650Yv2_manifest.csv.gz"
    ##   .. ..$ taxid                  : chr "9606"
    ##   .. ..$ technology             : chr "oligonucleotide beads"
    ##   .. ..$ title                  : chr "Illumina HumanHap650Yv2 Genotyping BeadChip (HumanHap650Y_v2-0_GenotypingBC_11237679_A)"

@dataTable contiene toda la informacion sobre la tabla de variantes del
chip. Dentro hay dos objetos @columns (que contiene información sobre
las columnas del chip) y @table que contiene los datos de las variantes
del chip:

``` r
tabla_datos_gpl<-gpl@dataTable
table_gpl<-(tabla_datos_gpl@table)
head(table_gpl)
```

    ##                                          ID        Name IlmnStrand SNP_ID
    ## 1 MitoA10045G-gi13273284_B_R_IFB1141652022: MitoA10045G        Bot  [T/C]
    ## 2 MitoA10551G-gi13273286_T_F_IFB1141639111: MitoA10551G        Top  [A/G]
    ## 3 MitoA11252G-gi13273288_B_R_IFB1141661584: MitoA11252G        Bot  [T/C]
    ## 4 MitoA11468G-gi13273289_B_R_IFB1141719645: MitoA11468G        Bot  [T/C]
    ## 5 MitoA11813G-gi13273292_T_F_IFB1141667833: MitoA11813G        Top  [A/G]
    ## 6 MitoA12309G-gi13273294_T_F_IFB1141664161: MitoA12309G        Top  [A/G]
    ##   AddressA_ID AlleleA_ProbeSeq AddressB_ID AlleleB_ProbeSeq Chr MapInfo Ploidy
    ## 1   904860368               NA           0               NA   M   10045   1000
    ## 2   905860131               NA           0               NA   M   10551   1000
    ## 3   903990168               NA           0               NA   M   11252   1000
    ## 4   901990681               NA           0               NA   M   11468   1000
    ## 5   904780630               NA           0               NA   M   11813   1000
    ## 6   906590470               NA           0               NA   M   12309   1000
    ##        Species CustomerStrand IlmnStrand IllumicodeSeq TopGenomicSeq
    ## 1 Homo sapiens            Top        Bot            NA            NA
    ## 2 Homo sapiens            Top        Top            NA            NA
    ## 3 Homo sapiens            Top        Bot            NA            NA
    ## 4 Homo sapiens            Top        Bot            NA            NA
    ## 5 Homo sapiens            Top        Top            NA            NA
    ## 6 Homo sapiens            Top        Top            NA            NA

El número de filas de este df es el numero de variantes del array:

``` r
nrow(table_gpl)
```

    ## [1] 662841

Tambien podemos verlo en el objeto @header:

``` r
gpl@header$data_row_count
```

    ## [1] "662841"

Comprobamos el número de individuos (muestras) del experimento:

``` r
length(gpl@header$sample_id)
```

    ## [1] 1215

Exploramos ahora una muestra para ver las variables que contiene el
experimento:

``` r
gsm <- getGEO("GSM895465")
gsm@header$characteristics_ch1
```

    ## [1] "qc status: Failed QC" "cell type: Blood"     "disease status: AD"

Esta muestra era de sangre y el individuo tenía alzheimer

Vemos además la tabla con los datos que tiene cada muestra:

``` r
head(Table(gsm))
```

    ##                                      ID_REF Allele1-Top Allele2-Top  Score
    ## 1 MitoA10045G-gi13273284_B_R_IFB1141652022:           A           A 0.3070
    ## 2 MitoA10551G-gi13273286_T_F_IFB1141639111:           A           A 0.3113
    ## 3 MitoA11252G-gi13273288_B_R_IFB1141661584:           A           A 0.7245
    ## 4 MitoA11468G-gi13273289_B_R_IFB1141719645:           A           A 0.3510
    ## 5 MitoA11813G-gi13273292_T_F_IFB1141667833:           A           A 0.5721
    ## 6 MitoA12309G-gi13273294_T_F_IFB1141664161:           A           A 0.6480
    ##   Theta     R B Allele Freq  VALUE
    ## 1 0.000 6.678             0 0.8744
    ## 2 0.007 7.450             0 1.1457
    ## 3 0.011 6.636             0 0.8919
    ## 4 0.003 9.025             0 0.4972
    ## 5 0.026 7.897             0 1.0520
    ## 6 0.007 7.181             0 1.0365
