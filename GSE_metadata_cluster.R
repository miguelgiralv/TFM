
# Comprobamos si tenemos GEOquery
sink("output_R.txt")

library(GEOquery)

setwd("/home/mgiralv/tfm/gse")

gse <- getGEO(filename = "/home/mgiralv/tfm/gse/files/GSE33528_family.soft.gz")

save(gse, file = "gse")

metadata_gse<-Meta(gse)
print(metadata_gse)

print("terminado!")

sink()



