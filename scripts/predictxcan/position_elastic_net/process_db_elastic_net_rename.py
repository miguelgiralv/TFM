import sqlite3
import pandas as pd
import os

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# primero vemos las tablas que hay en los archivos .db de predictdb
def list_tables(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    conn.close()
    return tables
list_tables(f"{path}/data/predictXcan/mashr/mashr_Whole_Blood.db")
# vemos que hay dos tablas, weights and extra
# hacemos una funcion para cargar las tablas como df de pandas
def fetch_data_as_df(db_path, table_name):
    conn = sqlite3.connect(db_path)
    query = f"SELECT * FROM {table_name}"
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

# cargamos la tabla de adipose subcutaneous 
weight_adipose_sub = fetch_data_as_df(f"{path}/data/predictXcan/mashr_models/mashr_Adipose_Subcutaneous.db", "weights")
weight_adipose_sub.head()
# vemos que hay un snp por fila, por lo que el modelo predictivo para este tejido será:
len(weight_adipose_sub)
#ahora hacemos esto para todas los tejidos del directorio y crearemos una tabla que resuma el numero de snps:
fetch_data_as_df(f"{path}/data/predictXcan/mashr_models/mashr_Adipose_Subcutaneous.db", "weights")
db_files = [f for f in os.listdir(f"{path}/data/predictXcan/mashr_models/") if f.endswith('.db')]
db_path=f"{path}/data/predictXcan/mashr_models"
tissues = []
n_snps = []
SNPS_tejidos=pd.DataFrame({"Tissue":[], "SNPs":[]})
for db_file in db_files:
    db_file_path_full =f"{db_path}/{db_file}"  
    db_df = fetch_data_as_df(db_file_path_full, "weights")
    n_snps.append(len(db_df))   
    table_name = os.path.splitext(db_file)[0]
    table_name_corrected = table_name[3:]
    tissues.append(table_name_corrected)  
    
SNPS_tejidos["Tissue"]=tissues
SNPS_tejidos["SNPs"]=n_snps


#ahora crearemos una tabla que contenga solo las posiciones de los snps de cada una de las tablas df
counter = 0
db_dataframes = {}

for db_file in db_files:
    counter += 1
    db_file_path_full = f"{db_path}/{db_file}"  
    db_dataframes[counter] = fetch_data_as_df(db_file_path_full, "weights")
print(counter)

total_df = pd.DataFrame({"varID": []})
for i in range(1, counter + 1):
    df_db = db_dataframes[i]
    rsid_column = df_db[['varID']]
    total_df = pd.concat([total_df, rsid_column]).drop_duplicates()
total_df

# Extraemos el cromosoma, la posicion y los valores de nucleotidosde la columna varID
total_df['CHROM'] = total_df['varID'].apply(lambda x: x.split('_')[0])
total_df['BEG'] = total_df['varID'].apply(lambda x: int(x.split('_')[1]))-1
total_df['END'] = total_df['varID'].apply(lambda x: int(x.split('_')[1]))

# Creating a new DataFrame with the extracted values
final_df = total_df[['CHROM', 'BEG', 'END']]

final_df.to_csv(f"{path}/results/imputado/extracted/all_positions_mashr.txt", sep="\t", header=False, index=False)

#con este archivo haremos un liftover a hg19 (filtrarvcf_position) y ahora le quitamos el chr para poder filtrar con bcftools:
input_file = f"{path}/results/imputado/extracted/all_positions_mashr_hg19.txt" 
output_file = f"{path}/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt" 

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        parts = line.strip().split('\t')  
        chromosome = parts[0].replace('chr', '') 
        start = parts[1]
        end = parts[2]
        outfile.write(f"{chromosome}\t{start}\t{end}\n") 

# y luego filtraremos los vcf con el  (filtrarvcf_position)

# Añadimos ahora esta información a nuestro anterior df  y luego lo guardamos todo
SNPs_total=len(total_df)
total_observado = {'Tissue': "Total", 'SNPs': SNPs_total}
total_observado_df = pd.DataFrame([total_observado])
SNPS_tejidos=pd.concat([SNPS_tejidos,total_observado_df])
SNPS_tejidos.to_csv(f"{path}/data/predictXcan/mashr_models/SNPs_tejidos_mashr.csv", index=False)

 
