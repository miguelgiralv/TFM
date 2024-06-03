import pandas as pd

#######################################################
path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# cargamos todos los datos (números de casos )
# Los números de casos antes del análisis de calidad en base a metadatos
with open(f"{path}/results/plink_data/binary/processed/GSE33528_sex.fam", "r") as file:
    row_count_qc_antes = len(file.readlines())    
#y despues del análisis de calidad en base a metadatos:
with open(f"{path}/results/plink_data/binary/processed/GSE33528_sex_qc.fam", "r") as file:
    row_count_qc_despues = len(file.readlines())
# y en base al call rate de individuos
with open(f"{path}/results/plink_data/binary/processed/GSE33528_qc_cr_in.fam", "r") as file:
    row_count_qc_in_despues = len(file.readlines())

# Resumimos toda la información:
print(" Antes de descartar los individuos con muestras de baja calidad según los metadatos, teníamos", row_count_qc_antes, "individuos.","\n",
"Tras el análisis de calidad inicial, nos quedamos con", row_count_qc_despues, "individuos.","\n",
"Hemos descartado", row_count_qc_antes-row_count_qc_despues, "individuos en el análisis de calidad de los metadatos.","\n","\n",
"Tras el análisis de call-rate, nos quedamos con", row_count_qc_in_despues, "individuos.","\n",
"Hemos descartado", row_count_qc_despues-row_count_qc_in_despues, "individuos.","\n","\n",
"En total se han descartado", row_count_qc_antes-row_count_qc_in_despues,"individuos")

# y guardamos todo en una tabla: 
qc_resumen=pd.DataFrame({
    "Inicial":[row_count_qc_antes],
    "QC metadatos": [row_count_qc_despues],
    "QC call rate": [row_count_qc_in_despues]     
})

qc_resumen.to_csv(f"{path}/results/metadata_gsm/qc_resumen.csv", index=False)

###########################################################
# Ahora haremos una análisis del número de variantes
# Primero veremos el número de SNPs que teníamos antes de mapear de hg18 a hg38
with open(f"{path}/results/plink_data/classic/GSE33528.map", "r") as file:
    row_count_var_map_antes = len(file.readlines())   
# y ahora vemos el numero de SNPs tras el mapeo:
with open(f"{path}/results/plink_data/classic/processed/output_map.map", "r") as file:
    row_count_var_map_despues = len(file.readlines())   
#Vemos el número de variantes que nos quedan despues de filtrar el call-rate:
with open(f"{path}/results/plink_data/binary/processed/GSE33528_qc_cr.bim", "r") as file:
    row_count_var_cr_despues = len(file.readlines())   

# mostramos toda la información:
print("Antes de mapear el genoma de referencia teníamos", row_count_var_map_antes, "variantes.","\n",
"Tras el mapeo a hg19, nos quedamos con", row_count_var_map_despues, "variantes.","\n",
"Hemos descartado", row_count_var_map_antes-row_count_var_map_despues, "variantes.","\n","\n",
"Tras el análisis de calidad en base al call-rate, nos quedamos con", row_count_var_cr_despues, "variantes.","\n",
"Hemos eliminado", row_count_var_map_despues-row_count_var_cr_despues, "variantes.","\n","\n",
"En total se han descartado", row_count_var_map_antes-row_count_var_cr_despues,"variantes")

# ahora añadimos el numero de variantes que quedan despues de imputar:
imputado=pd.read_csv(f"{path}/results/imputado/extracted/merged.info", compression='gzip')

SNPs=sum(imputado["Genotyped"]=="Genotyped")

SNPs_total=len(imputado)

# y la resumimos en una tabla que guardamos:
variantes_resumen=pd.DataFrame({
    "Inicial":[row_count_var_map_antes],
    "Mapeo": [row_count_var_map_despues],
    "Call-rate": [row_count_var_cr_despues],
    "Imputado_originales": [SNPs],
    "Imputado":SNPs_total         
})

variantes_resumen.to_csv(f"{path}/results/metadata_gsm/variantes_resumen.csv", index=False)
