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
list_tables(f"{path}/data/predictXcan/elastic_net_models/en_Adipose_Subcutaneous.db")
# vemos que hay dos tablas, weights and extra
# hacemos una funcion para cargar las tablas como df de pandas
def fetch_data_as_df(db_path, table_name):
    conn = sqlite3.connect(db_path)
    query = f"SELECT * FROM {table_name}"
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

# cargamos la tabla de adipose subcutaneous 
weight_adipose_sub = fetch_data_as_df(f"{path}/data/predictXcan/elastic_net_models/en_Adipose_Subcutaneous.db", "weights")
weight_adipose_sub.head()
# vemos que hay un snp por fila, por lo que el modelo predictivo para este tejido será:
len(weight_adipose_sub)
#ahora hacemos esto para todas los tejidos del directorio y crearemos una tabla que resuma el numero de snps:
fetch_data_as_df(f"{path}/data/predictXcan/elastic_net_models/en_Adipose_Subcutaneous.db", "weights")
db_files = [f for f in os.listdir(f"{path}/data/predictXcan/elastic_net_models") if f.endswith('.db')]
db_path=f"{path}/data/predictXcan/elastic_net_models"
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


#ahora crearemos una tabla que contenga solo los rsids de cada una de las tablas df
counter = 0
db_dataframes = {}

for db_file in db_files:
    counter += 1
    db_file_path_full = f"{db_path}/{db_file}"  
    db_dataframes[counter] = fetch_data_as_df(db_file_path_full, "weights")
print(counter)

total_df = pd.DataFrame({"rsid": []})
for i in range(1, counter + 1):
    df_db = db_dataframes[i]
    rsid_column = df_db[['rsid']]
    total_df = pd.concat([total_df, rsid_column]).drop_duplicates()
total_df

total_df.to_csv(f"{path}/results/predictxcan/all_rsids.txt", sep="\t", header=False, index=False)

# Añadimos ahora esta información a nuestro anterior df  y luego lo guardamos todo
SNPs_total=len(total_df)
total_observado = {'Tissue': "Total", 'SNPs': SNPs_total}
total_observado_df = pd.DataFrame([total_observado])
SNPS_tejidos=pd.concat([SNPS_tejidos,total_observado_df])
SNPS_tejidos.to_csv(f"{path}/data/predictXcan/elastic_net_models/SNPs_tejidos.csv", index=False)

 
# con esta lista filtraremos luego los vcf de los cromosomas



   

