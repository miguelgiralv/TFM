import pandas as pd
import os


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

# cogemos los SNPs de la bd de predictDB con la tabla weights

weights=pd.read_csv(f"{path}/data/predictXcan/models/gtex_v8_en/exported_weights_db.csv")

# Cogemos y procesamos la columna varID de la tabla weights, los procesamos: chr10_295356_T_C_b38  -> 10:295356 
# y los guardamos en un txt de aquí los seleccionaremos en el vcf con bcftools para eliminar los que no están en 

def process_snp_id(snp_id):
    parts = snp_id.split('_')  # Split the SNP ID by underscores
    chromosome = parts[0].replace('chr', '')  # Remove 'chr' prefix from chromosome
    position = parts[1]  # Position remains unchanged
    return f"{chromosome}:{position}"  # Format the output

weights["varID"] = weights["varID"].apply(process_snp_id)
len(weights)
# tenemos 206238  snps en predictdb
weights["varID"].to_csv(f"{path}/results/imputado/extracted/snps_predictdb.txt", sep="\t", index=False, header=False)




