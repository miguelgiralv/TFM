
import pandas as pd
# Cargamos la tabla de fenotipos que hemos creado (en metadata_table_generate) con los metadatos 
path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

metadata = pd.read_csv(f"{path}/results/metadata_gsm/merged_metadata.txt", sep="\t")
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

# Guardamos la informacion en una tabla para cargarlo en el informe 
tabla_resumen_casos = pd.DataFrame({
    "Total": [n_total],
    "Control": [n_control],
    "Casos": [n_AD],
    "Desconocido":[n_unknown]
})
tabla_resumen_casos.to_csv(f"{path}/results/metadata_gsm/tabla_resumen_casos.csv",index=False)

print("En un inicio tenemos",n_total, "individuos,","de los cuales" "\033[92m",n_AD,"\033[0m", "tienen Alzheimer,","\033[33m", n_control, "\033[0m" "son individuos control y", 
"\033[91m", n_unknown, "\033[0m","son individuos con fenotipo desconocido.\n")

tabla_sexos = pd.DataFrame({
    "Total": [n_total],
    "Hombre": [n_male],
    "Mujeres": [n_female],
    "Desconocido":[n_unknown]
})

tabla_sexos.to_csv(f"{path}/results/metadata_gsm/tabla_sexos.csv",index=False)

print("En un inicio tenemos,", n_male,"hombres y", n_female, " mujeres.\n Por otro lado","\033[92m", n_qc_passed, "\033[0m","de las muestras pasaron el control de calidad y",
"\033[91m", n_qc_failed, "\033[0m","de las no lo pasaron")

# tras realizar todos los análisis de calidad, hemos descartado varias muestras, así que ahora veremos las características fenotípicas de estas muestras procesadas
# Primero obtenemos el df que contiene a los individuos con los que nos hemos quedado despues de todos los análisis de calidad (en el archivo fam):

print(f"{path}/results/plink_data/binary/processed/GSE33528_qc_cr_in.fam")
fam_df=pd.read_csv(f"{path}/results/plink_data/binary/processed/GSE33528_qc_cr_in.fam",
names=["fam_id", "Individual ID", "1","2", "3","4"], sep=r"\s+")

print(fam_df.head(f"{path}/results/metadata_gsmprocessed/GSE33528_qc_cr_in.fam"))

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

# generamos las tablas:
tabla_resumen_casos_qc = pd.DataFrame({
    "Total": [n_total_qc],
    "Control": [n_control_qc],
    "Casos": [n_AD_qc],
    "Desconocido":[n_unknown_qc]
})
tabla_resumen_casos_qc.to_csv(f"{path}/results/metadata_gsm/tabla_resumen_casos_qc.csv",index=False)

tabla_sexos_qc = pd.DataFrame({
    "Total": [n_total],
    "Hombres": [n_male_qc],
    "Mujeres": [n_female_qc]
})
tabla_sexos_qc.to_csv(f"{path}/results/metadata_gsm/tabla_sexos_qc.csv",index=False)

# y hacemos un recuento:
print(" Tras el análisis de calidad, tenemos", n_total_qc, "individuos,","de los cuales" "\033[92m",n_AD_qc,"\033[0m", "tienen Alzheimer,","\033[33m", n_control_qc, "\033[0m" "son individuos control y", 
"\033[91m", n_unknown_qc, "\033[0m","son individuos con fenotipo desconocido.\n","De nuestros individuos,", n_male_qc,"son hombres y", n_female_qc, "son mujeres")



