import pandas as pd
import GEOparse
import re
import os

############################################
# Generamos una tabla con los metadatos:
def GSM_metadata (gsm_inicial, gsm_final):
    dfs = []
    for i in range(gsm_inicial, gsm_final+1):
        gsm_accession="GSM" + str(i)
        gsm = GEOparse.get_GEO(gsm_accession, destdir= "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/GSMs")
        # Extraemos si pasó el control de calidad primero
        characteristics = gsm.metadata["characteristics_ch1"]
        if characteristics[0] == "qc status: Passed QC":
            qc = "Passed"
        elif characteristics[0] == "qc status: Failed QC":
            qc = "Failed"
        else:
            qc = "unknown"    
        # Extraemos el tipo celular de los metadatos
        match = re.search(r'cell type:\s*(\w+)', characteristics[1], re.IGNORECASE)
        if match:
            cell_type = match.group(1)
        else:
            cell_type = "unknown"
        # Extraemos el phenotipo de la enfermedad los metadatos
        if characteristics[2] == "disease status: Unaffected":
            phenotype = "control"
        elif characteristics[2] == "disease status: AD":
            phenotype = "AD"
        else:
            phenotype = "unknown"       
        new_row_df = pd.DataFrame({
            "Individual ID": gsm_accession,
            "QC": qc,
            "Cell type": cell_type,
            "Disease status": phenotype,
            "Individual ID": gsm_accession
        }, index=[0])
        # actualizamos la lista que almacena las filas del df
        dfs.append(new_row_df)
    # generamos el df final
    result_df = pd.concat(dfs, ignore_index=True)
    print(result_df)
    return(result_df)

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
# test
df0=GSM_metadata(894502, 894503)
#generamos varios dfs (algunos gsm ids en medio no son del mismo experimento)
df1=GSM_metadata(894502, 894572).to_csv(f"{path}/results/metadata_gsm/metadatas/metadata_1.txt", index=False, header=True, sep="\t")
df2=GSM_metadata(894589, 894801).to_csv(f"{path}/results/metadata_gsm/metadatas/metadata_2.txt", index=False, header=True, sep="\t")
df3=GSM_metadata(894825, 894895).to_csv(f"{path}/results/metadata_gsm/metadatas/metadata_3.txt", index=False, header=True, sep="\t")
df4=GSM_metadata(894897, 894967).to_csv(f"{path}/results/metadata_gsm/metadatas/metadata_4.txt", index=False, header=True, sep="\t")
df5=GSM_metadata(894980, 895768).to_csv(f"{path}/results/metadata_gsm/metadatas/metadata_5.txt", index=False, header=True, sep="\t")

# unimos los dfs en uno único que contiene todos los metadatos:
merged_metadata = pd.DataFrame()
metadata_dir = f"{path}/results/metadata_gsm/metadatas/"
for metadata_file in os.listdir(metadata_dir):
    file_path = os.path.join(metadata_dir, metadata_file)
    df = pd.read_csv(file_path, delimiter="\t") 
    merged_metadata = pd.concat([merged_metadata, df], ignore_index=True)

#comprobamos la tabla de metadatos:
merged_metadata.head
# vemos si hay columnas redundantes (con un solo valor, tal como cell type=blood) y las eliminamos:
for column_name in merged_metadata.columns:
    if merged_metadata[column_name].nunique() <2:
        merged_metadata=merged_metadata.drop(columns=column_name)
# hemos eliminado cell type porque solo tiene un valor
merged_metadata.head
# a continuación añadimos el sexo a nuestra tabla de metadatos, para ello usaremos la tabla checksex que generamos
# anteriormente al corregir el sexo con la funcion check-sex:
check_sex=pd.read_csv(f"{path}/results/metadata_gsm/GSE33528_sex_check_results.sexcheck", sep="\s+")
check_sex.head()
# armonizamos nombres en la columna de individuos de checksex 
check_sex.rename(columns={'IID':'Individual ID'}, inplace=True)
# y ahora fusionamos los dfs
merged_df = pd.merge(merged_metadata, check_sex, on='Individual ID', how='inner')
# nos quedamos con las columnas de interes
merged_metadata = merged_df[['Individual ID', 'QC', 'Disease status', 'PEDSEX']]
# cambiamos nombre de la columna sexo:
merged_metadata.rename(columns={'PEDSEX':'Sex'}, inplace=True)
# sustituimos los numeros por male y female
merged_metadata['Sex'] = merged_metadata['Sex'].replace({1: 'male', 2: 'female'})

# guardamos finalmente la tabla con metadatos:
merged_metadata.to_csv(f"{path}/results/metadata_gsm/merged_metadata.txt", index=False, header=True, sep="\t")
