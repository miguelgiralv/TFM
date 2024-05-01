
import pandas as pd
# Cargamos la tabla de fenotipos que hemos creado (en metadata_table_generate) con los metadatos 
metadata = pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", sep="\t")
# La vemos
metadata

# Observamos las características fenotípicas antes de hacer el análisis de calidad 
n_total=metadata.shape[0]
n_qc_passed=sum(metadata["QC"]=="Passed")
n_qc_failed=sum(metadata["QC"]=="Failed")
n_AD=sum(metadata["Disease status"]=="AD")
n_control=sum(metadata["Disease status"]=="control")
n_unknown=sum(metadata["Disease status"]=="unknown")
n_male=sum(metadata["Sex"]=="male")
n_female=sum(metadata["Sex"]=="female")
# y hacemos un recuento:
print(" En un inicio tenemos",n_total, "individuos,","de los cuales" "\033[92m",n_AD,"\033[0m", "tienen Alzheimer,","\033[33m", n_control, "\033[0m" "son individuos control y", 
"\033[91m", n_unknown, "\033[0m","son individuos con fenotipo desconocido.\n","De estos,", n_male,"son hombres y", n_female, "son mujeres.\n Por otro lado","\033[92m", n_qc_passed, "\033[0m","de las muestras pasaron el control de calidad y",
"\033[91m", n_qc_failed, "\033[0m","de las no lo pasaron")

# tras realizar todos los análisis de calidad, hemos descartado varias muestras, así que ahora veremos las características fenotípicas de estas muestras procesadas
# Primero obtenemos el df que contiene a los individuos con los que nos hemos quedado despues de todos los análisis de calidad (en el archivo fam):
fam_df=pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.fam",
names=["fam_id", "Individual ID", "1","2", "3","4"], sep="\s+")
# cogemos los individuos que han pasado los análisis de calidad con set
individual_ids = set(fam_df['Individual ID'])
# y filtramos en nuestra tabla de metadatos:
metadata_qc = metadata[metadata['Individual ID'].isin(individual_ids)]

# y ahora resumimos los fenotipos resultantes:
n_total_qc=metadata_qc.shape[0]
n_qc_passed_qc=sum(metadata_qc["QC"]=="Passed")
n_qc_failed_qc=sum(metadata_qc["QC"]=="Failed")
n_AD_qc=sum(metadata_qc["Disease status"]=="AD")
n_control_qc=sum(metadata_qc["Disease status"]=="control")
n_unknown_qc=sum(metadata_qc["Disease status"]=="unknown")
n_male_qc=sum(metadata_qc["Sex"]=="male")
n_female_qc=sum(metadata_qc["Sex"]=="female")
# y hacemos un recuento:
print(" Tras el análisis de calidad, tenemos", n_total_qc, "individuos,","de los cuales" "\033[92m",n_AD_qc,"\033[0m", "tienen Alzheimer,","\033[33m", n_control_qc, "\033[0m" "son individuos control y", 
"\033[91m", n_unknown_qc, "\033[0m","son individuos con fenotipo desconocido.\n","De nuestros individuos,", n_male_qc,"son hombres y", n_female_qc, "son mujeres")