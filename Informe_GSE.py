#######################################################################
# Control del numero de muestras y variantes que vamos teniendo con el análisis de calidad

#Vemos el número de muestras que teniamos antes del qc:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528.fam", "r") as file:
    row_count_qc_antes = len(file.readlines())    
#y despues del qc
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc.fam", "r") as file:
    row_count_qc_despues = len(file.readlines())    
# y despues de filtrar individuos con call-rate menor:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr.fam", "r") as file:
    row_count_qc_cr_despues = len(file.readlines())   

print("Antes del análisis de calidad inicial, teníamos", row_count_qc_antes, "individuos.","\n",
"Tras el análisis de calidad inicial, nos quedamos con", row_count_qc_despues, "individuos.","\n",
"Hemos eliminado", row_count_qc_antes-row_count_qc_despues, "individuos en el análisis de calidad de los metadatos.","\n",
"Tras el filtrado por call rate, nos quedamos", row_count_qc_cr_despues, "individuos.","\n",
"En este paso de filtrado, hemos eliminado", row_count_qc_despues-row_count_qc_cr_despues, "individuos")

#Vemos el número de variantes que teniamos antes del call-rate:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528.bim", "r") as file:
    row_count_cr_antes = len(file.readlines())    
#y despues de filtrar el call-rate:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr.bim", "r") as file:
    row_count_cr_despues = len(file.readlines())    
print("Inicialmente teníamos", row_count_cr_antes, "variantes.""\n",
"Tras filtrar variantes con call rate >0.95, nos quedamos", row_count_cr_despues, "variantes","\n",
"Hemos descartado", row_count_cr_antes-row_count_cr_despues, "variantes")