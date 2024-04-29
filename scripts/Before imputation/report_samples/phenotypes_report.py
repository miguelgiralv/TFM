import pandas as pd
import subprocess

######################################################################
# Cargamos la tabla de fenotipos que hemos creado con los metadatos 
merged_metadata = pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", sep="\t")
# La vemos
merged_metadata
# Vemos el número de individuos que tenemos primero en el archivo de plink binario
plink_command=["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink",
"--bfile",
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/GSE33528",
"--freq"]  
subprocess.run(plink_command)

# comprobamos la consistencia de nuestro metodo de asignación de sexo con el comando --check-sex de plink:
plink_command = [
    "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink",
    "--bfile",
    "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/GSE33528",
    "--check-sex",
    "--out",
    "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/GSE33528_sex_check_results"
]
subprocess.run(plink_command)
# abrimos el archivo generado y vemos la consistencia del sexo
sex_consistency=pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/GSE33528_sex_check_results.sexcheck", sep="\s+")
sex_ok = sex_consistency[sex_consistency['STATUS'] == 'OK'].shape[0]
sex_problem = sex_consistency[sex_consistency['STATUS'] == 'PROBLEM'].shape[0]
print(sex_ok,"individuos concuerdan con nuestro método,", sex_problem,"individuos no concuerdan con la predicción de plink")
# para arreglar estos 43 individuos, vamos a ejecutar --check-sex pero con parámetros un poco mas laxos:


# Cargamos la tabla de fenotipos que hemos creado con los metadatos 
merged_metadata = pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", sep="\t")
# Comprobamos que las dimensiones coinciden con el número de muestras
merged_metadata
# Vemos los distintos fenotipos en base a los metadatos:
control_count = merged_metadata[merged_metadata['Disease status'] == 'control'].shape[0]
AD_count = merged_metadata[merged_metadata['Disease status'] == 'AD'].shape[0]
unknown_count = merged_metadata[merged_metadata['Disease status'] == 'unknown'].shape[0]
qc_failed = merged_metadata[merged_metadata['QC'] == 'Failed'].shape[0]
qc_pass = merged_metadata[merged_metadata['QC'] == 'Passed'].shape[0]
# Y hacemos un recuento 
print("Tenemos", "\033[92m", AD_count,"\033[0m", "individuos con Alzheimer,","\033[33m", control_count, "\033[0m" "individuos control y", "\033[91m", unknown_count, 
"\033[0m","individuos con fenotipo desconocido.\n", "Por otro lado","\033[92m", qc_pass, "\033[0m","de las muestras pasaron el control de calidad y",
"\033[91m", qc_failed, "\033[0m","muestras fallaron el control de calidad")
 
plink_command=["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink",
"--bfile",
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/GSE33528",
"--freq"]  
subprocess.run(plink_command)
plink --file your_data --check-sex --out sex_check_results







