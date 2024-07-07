import pandas as pd
import os

# con este scripts sacaremos el numero de snps que se tienen despues de imputar en el michigan. 
# el resutlado de este script se usa en el informe QC_report
path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

info_dir = f"{path}/results/imputado/extracted"

info_dfs = []

for chrom in range(1, 23):
    filepath = os.path.join(info_dir, f"chr{chrom}.info.gz")  
    info_df = pd.read_csv(filepath, sep="\t",  on_bad_lines='skip')
    # almacenamos los valores numericos
    # lo almacenamos en una lista de dfs 
    info_dfs.append(info_df)

# luego metemos el cromosoma x
filepath_X = os.path.join(info_dir, "chrX.info.gz")
info_df_X = pd.read_csv(filepath_X, sep="\t",  on_bad_lines='skip')
info_dfs.append(info_df_X)

# ahora concatenamos nuestra lista de dfs para obtener un unico df
all_info_data = pd.concat(info_dfs, ignore_index=True)

# guardamos el df
all_info_data.to_csv(f"{path}/results/imputado/extracted/merged.info", index=False, compression='gzip')

# vemos cuantos SNPs son genotyped
sum(all_info_data["Genotyped"]=="Genotyped")
# 641576 SNPs
sum(all_info_data["Genotyped"]=="Imputed")

########### hacemos lo mismo con los datos de imputacion con el parametro 0,1
info_dir_2 = f"{path}/results/imputado_2"

info_dfs_2 = []

for chrom in range(1, 23):
    filepath = os.path.join(info_dir_2, f"chr{chrom}.info.gz")  
    info_df_2 = pd.read_csv(filepath, sep="\t",  on_bad_lines='skip')
    # almacenamos los valores numericos
    # lo almacenamos en una lista de dfs 
    info_dfs_2.append(info_df_2)

# luego metemos el cromosoma x
filepath_X_2 = os.path.join(info_dir_2, "chrX.info.gz")
info_df_X_2 = pd.read_csv(filepath_X_2, sep="\t",  on_bad_lines='skip')
info_dfs_2.append(info_df_X_2)

# ahora concatenamos nuestra lista de dfs para obtener un unico df
all_info_data_2 = pd.concat(info_dfs_2, ignore_index=True)


# guardamos el df
all_info_data_2.to_csv(f"{path}/results/imputado/extracted/merged_2.info", index=False, compression='gzip')

# vemos cuantos SNPs son genotyped
sum(all_info_data_2["Genotyped"]=="Genotyped")
# 641576 SNPs
sum(all_info_data_2["Genotyped"]=="Imputed")



