import pandas as pd


#######################################################
# Primero hacemos un recuento de individuos a lo largo de todo el análisis:
# Eliminaremos primero los individuos con qc=failed en sus metadatos
merged_metadata = pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", sep="\t")
qc_failed=merged_metadata[merged_metadata["QC"]=="Failed"]
#creamos un df con la lista de ids de muestras a eliminar y su family ID, que es siempre el mismo ()
qc_failed_ids = qc_failed[["Individual ID"]]
qc_failed_ids['Family_id'] = "GSE33528"
# cambiamos el orden para que family_id sea la primera columna
qc_failed_ids = qc_failed_ids[['Family_id', 'Individual ID']]
# lo almacenamos en un txt
qc_failed_ids.to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt", sep="\t", index=False, header=False, lineterminator="\n")
# ahora ejecutamos plink para eliminar las muestras que fallaron el test de calidad y que hemos almacenado en qc_failed_ids.txt
import subprocess
# el comando es:
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/GSE33528", 
"--remove", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt","--allow-extra-chr", 
"--make-bed", "--out", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc"]
# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)
# antes del qc:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528.fam", "r") as file:
    row_count_qc_antes = len(file.readlines())    
#y despues del qc:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc.fam", "r") as file:
    row_count_qc_despues = len(file.readlines())   

print("Antes de hacer descartar los individuos con muestras de baja calidad según los metadatos, teníamos", row_count_qc_antes, "individuos.","\n",
"Tras el análisis de calidad inicial, nos quedamos con", row_count_qc_despues, "individuos.","\n",
"Hemos eliminado", row_count_qc_antes-row_count_qc_despues, "individuos en el análisis de calidad de los metadatos.")

# Ahora seleccionamos ahora los individuos con buena calidad en base al call-rate de individuos (comando --mind)
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc", 
"--mind", "0.05","--allow-extra-chr","--out",
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_in",
"--make-bed"] 
# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_in.fam", "r") as file:
    row_count_qc_in_despues = len(file.readlines())   

# Y resumimos toda la información:

print("Antes de descartar los individuos con muestras de baja calidad según los metadatos, teníamos", row_count_qc_antes, "individuos.","\n",
"Tras el análisis de calidad inicial, nos quedamos con", row_count_qc_despues, "individuos.","\n",
"Hemos descartado", row_count_qc_antes-row_count_qc_despues, "individuos en el análisis de calidad de los metadatos.","\n","\n",
"Tras el análisis de call-rate, nos quedamos con", row_count_qc_in_despues, "individuos.","\n",
"Hemos descartado", row_count_qc_despues-row_count_qc_in_despues, "individuos.","\n","\n",
"En total se han descartado", row_count_qc_antes-row_count_qc_in_despues,"individuos")

###########################################################
# Ahora haremos una análisis del número de variantes
# Primero veremos el número de SNPs que teníamos antes de mapear de hg18 a hg38
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map", "r") as file:
    row_count_var_map_antes = len(file.readlines())   

# y ahora vemos el numero de SNPs tras el mapeo:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/processed/output_map.map", "r") as file:
    row_count_var_map_despues = len(file.readlines())   

# Seleccionamos ahora solo las variantes con buena calidad en base a su call-rate (call rate de 0.95) (comando --geno)
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_in", 
"--geno", "0.05","--allow-extra-chr",
"--out", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final",
"--make-bed"] 
# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)
#Vemos el número de variantes que nos quedan despues de filtrar el call-rate:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.bim", "r") as file:
    row_count_var_cr_despues = len(file.readlines())   

# y resumimos toda la información:
print("Antes de mapear el genoma de referencia teníamos", row_count_var_map_antes, "variantes.","\n",
"Tras el mapeo a hg38, nos quedamos con", row_count_var_map_despues, "variantes.","\n",
"Hemos descartado", row_count_var_map_antes-row_count_var_map_despues, "variantes.","\n","\n",
"Tras el análisis de calidad en base al call-rate, nos quedamos con", row_count_var_cr_despues, "variantes.","\n",
"Hemos eliminado", row_count_var_map_despues-row_count_var_cr_despues, "variantes.","\n","\n",
"En total se han descartado", row_count_var_map_antes-row_count_var_cr_despues,"variantes")