#map rsid to position
library(biomaRt)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/" # el path al repositorio
file_path="data/PRS/PGS000054_hmPOS_GRCh37_2.txt"


rsid_1 <- read.csv(paste0(path,file_path), header = TRUE, sep = "\t")

rsids<-rsid_1$rsID

class(rs)

ensembl <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org")


snp_data <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start", "allele"),
  filters = "snp_filter",
  values = rsids,
  mart = ensembl
)


### FUNCIONAAA

class(snp_data)
library(dplyr)
library(tidyr)

snp_data <- snp_data %>%
  separate(allele, into = c("ref", "alt"), sep = "/") %>%
  mutate(
    formatted = paste0(chr_name, ":", chrom_start, ":", ref, ":", alt)
  )



snp_data <- snp_data %>%
  mutate(
    chr_name = ifelse(grepl("^HSCHR", chr_name), 
                      str_extract(chr_name, "\\d+"), 
                      as.character(chr_name))
  )


# ajustar, algunos rsids tienen varios alelos de referencia, ver como estan en mi vcf

snp_data

nt.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start",                   
                      "chrom_strand","associated_gene",
                      "ensembl_gene_stable_id"),
                    filters="refsnp",
                    values=rsids,
                    mart=ensembl)


listAttributes(ensembl)
