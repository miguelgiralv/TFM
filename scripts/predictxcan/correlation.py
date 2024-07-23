import pandas as pd
import os
import numpy as np

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

files_list = os.listdir(tissue_paths)

predict_files = [file for file in files_list if file.endswith("predict.txt")]

data_frames = []

tissue=predict_files[1]
for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID"], axis=1)
    # quitamos el indice que esta ahora como headers
    tissue_name = tissue.replace("mashr_", "").replace("predict.txt", "")
    # Renombramos las columnas para que cada gen una tenga el nombre del tejido
    tissue_pd.columns = [f"{col}_{tissue_name}" if col != 'FID' else col for col in tissue_pd.columns]
    data_frames.append(tissue_pd)
    print(tissue,"cargado")

###### buenoo
tissue=predict_files[1]
tissue_pd=pd.DataFrame

data_frames = []

for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID"], axis=1)
    tissue_pd=tissue_pd.T
    # quitamos el indice que esta ahora como headers
    tissue_pd.columns = tissue_pd.iloc[0]
    tissue_pd = tissue_pd[1:]
    tissue_name = tissue.replace("mashr_", "").replace("predict.txt", "")
    # Renombramos las columnas para que cada gen una tenga el nombre del tejido
    tissue_pd.columns = [f"{col}_{tissue_name}" if col != 'IID' else col for col in tissue_pd.columns]
    data_frames.append(tissue_pd)
    
    print(tissue,"cargado")
  
merged_data = pd.concat(data_frames, axis=1, join='inner')
print(merged_data.head())

    filepath = os.path.join(info_dir, f"chr{chrom}.info.gz")  
    info_df = pd.read_csv(filepath, sep="\t",  on_bad_lines='skip')
    # almacenamos los valores numericos
    # lo almacenamos en una lista de dfs 
    info_dfs.append(info_df)


merged_data = df_tissue1.join(df_tissue2, lsuffix='_tissue1', rsuffix='_tissue2')


prs_file = pd.read_csv(f"{path}/data/PRS/PGS004146_hmPOS_GRCh37_2.txt", sep="\t")

data_frames = [pd.read_csv(file) for file in file_paths]

for tissue in range(1, 23):
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



import pandas as pd
import os

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test/vcf_1000G_hg37_mashr"
files_list = os.listdir(path)
predict_files = [file for file in files_list if file.endswith("predict.txt")]

data_frames = []
for tissue in predict_files:
    # Load the DataFrame
    tissue_pd = pd.read_csv(os.path.join(path, tissue), sep="\t", on_bad_lines='skip')
    
    print(f"Processing file: {tissue}")
    print(tissue_pd.head())  # Check if data is loaded correctly
    
    # Ensure 'IID' is correctly set as index
    tissue_pd.set_index('IID', inplace=True)

    # Drop 'FID' if it exists (it seems to be redundant now)
    if 'FID' in tissue_pd.columns:
        tissue_pd = tissue_pd.drop(['FID'], axis=1)
    
    # Extract tissue name
    tissue_name = tissue.replace("mashr_", "").replace("predict.txt", "")

    # Rename columns to reflect tissue name
    tissue_pd.columns = [f"{col}_{tissue_name}" for col in tissue_pd.columns]

    # Append to list of DataFrames
    data_frames.append(tissue_pd)
    print(tissue_pd.head())  # Check if data is loaded correctly
    print(f"{tissue} loaded")


# Merge all DataFrames on 'FID'
merged_data = pd.concat(data_frames, axis=1, join='inner')
print(merged_data.head())
