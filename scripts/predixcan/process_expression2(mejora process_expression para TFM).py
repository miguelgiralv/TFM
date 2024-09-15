import os
import pandas as pd

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
### APLICAR A ELASTIC NET

summary_files = [f for f in os.listdir(os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en")) if f.endswith('summary.txt')]

# creamos un diccionario vacio para almacenar los df filtrados
filtered_dataframes = {}

# hacemos un bucle for en todos los archivos summary
for summary_file in summary_files:
    # leemos cada archivo summary:
    predict_file = summary_file.replace('_summary.txt', '_predict.txt')
    
    summary_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en", summary_file)
    predict_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en", predict_file)
    
    df_summary = pd.read_csv(summary_file_path, sep="\t")
    
    # seleccionamos todos los genes que tengan el valor NA en la columna de n_snps_used
    genes_to_remove = df_summary.loc[df_summary['n_snps_used'].isna(), 'gene'].tolist()
    
    # Ahora cargamos el archivo predict con los valores de expresion de los genes
    df_predict = pd.read_csv(predict_file_path, sep="\t")
    
    # guardamos los genes que hay que conservar (las columnas del archivo predict) en base a los genes con valor NA, tambien quitamos 'FID', porque es redundante
    columns_to_keep = [col for col in df_predict.columns if col not in genes_to_remove and col != 'FID']
    
    # Finalmente filtramos el df del archivo predict y lo guardamos 
    filtered_df_predict = df_predict[columns_to_keep]
    
    path_save = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en/processed", f"processed_{predict_file}")
    filtered_df_predict.to_csv(path_save, sep="\t", index=False)
    
    print(f"Processed {predict_file}: Initial dimensions={df_predict.shape}, Filtered dimensions={filtered_df_predict.shape}")
    
    filtered_dataframes[summary_file] = filtered_df_predict

   
### APLICAR A MASHR
summary_files = [f for f in os.listdir(os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr")) if f.endswith('summary.txt')]

filtered_dataframes = {}

for summary_file in summary_files:
    predict_file = summary_file.replace('_summary.txt', '_predict.txt')
    
    summary_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr", summary_file)
    predict_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr", predict_file)
    
    df_summary = pd.read_csv(summary_file_path, sep="\t")
    
    genes_to_remove = df_summary.loc[df_summary['n_snps_used'].isna(), 'gene'].tolist()
    
    df_predict = pd.read_csv(predict_file_path, sep="\t")
    
    columns_to_keep = [col for col in df_predict.columns if col not in genes_to_remove and col != 'FID']
    
    filtered_df_predict = df_predict[columns_to_keep]
    
    path_save = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr/processed", f"processed_{predict_file}")
    
    filtered_df_predict.to_csv(path_save, sep="\t", index=False)
    
    print(f"Processed {predict_file}: Initial dimensions={df_predict.shape}, Filtered dimensions={filtered_df_predict.shape}")    
    filtered_dataframes[summary_file] = filtered_df_predict


### Ahora obtenemos el numero de genes por cada tejido para elastic net y mashr antes y despues de procesarse:

# antes de procesarse (cogemos el summary.txt)
def count_genes_in_files(custom_path, file_prefix, model_name):
    predict_files = [f for f in os.listdir(custom_path) if f.startswith(file_prefix) and f.endswith("_summary.txt")]
    gene_list = []    # lista que almacena todos los genes
    genes_df = {}

    for predict_file in predict_files:
        start_index_name = predict_file.find(file_prefix) + len(file_prefix)
        end_index_name = predict_file.find("_summary.txt")
        predict_file_name = predict_file[start_index_name:end_index_name]       
        predict_file_path = os.path.join(custom_path, predict_file)
        df_predict = pd.read_csv(predict_file_path, sep="\t")
        n_genes = len(df_predict) - 1  # el encabezado
        gene_sublist=df_predict["gene"].tolist()
        gene_list= gene_list + gene_sublist        
        genes_df[predict_file_name] = n_genes
        print(predict_file, "finished")

    genes_df_final = pd.DataFrame(genes_df.items(), columns=['Tejido', model_name])
    unique_gene_list = list(set(gene_list)) # quitamos los genes repetidos de nuestra lista de genes para tener el total de genes que hay en todos los tejidos
    len_unique_gene_list=len(unique_gene_list) # numero de genes unicos en total usados en todos los tejidos
    genes_df_final.loc[len(genes_df_final)] = ['Total', len_unique_gene_list]
    
    return genes_df_final

def count_genes_in_files_processed(custom_path, file_prefix, model_name):
    gene_list = []
    predict_files = [f for f in os.listdir(custom_path) if f.startswith(file_prefix) and f.endswith("_predict.txt")]
    genes_df = {}
    for predict_file in predict_files:
        start_index_name = predict_file.find(file_prefix) + len(file_prefix)
        end_index_name = predict_file.find("_predict.txt")
        predict_file_name = predict_file[start_index_name:end_index_name]        
        predict_file_path = os.path.join(custom_path, predict_file)
        df_predict = pd.read_csv(predict_file_path, sep="\t")
        gene_sublist=(df_predict.columns.tolist())
        gene_list= gene_list + gene_sublist
        n_genes = len(df_predict.columns) - 1  # asumiendo que no hemos procesado el archivo aun (columna IID)
        genes_df[predict_file_name] = n_genes
        print(predict_file, "finished")

    genes_df_final = pd.DataFrame(genes_df.items(), columns=['Tejido', model_name])
    unique_gene_list = list(set(gene_list)) # quitamos los genes repetidos de nuestra lista de genes para tener el total de genes que hay en todos los tejidos
    len_unique_gene_list=len(unique_gene_list) # numero de genes unicos en total usados en todos los tejidos
    genes_df_final.loc[len(genes_df_final)] = ['Total', len_unique_gene_list]
    
    return genes_df_final

mashr_n_genes_before=count_genes_in_files("results/predictXcan/test/vcf_1000G_hg37_mashr", "mashr_", "Genes en mashr")
mashr_n_genes_processed=count_genes_in_files_processed("results/predictXcan/test/vcf_1000G_hg37_mashr/processed", "processed_mashr_", "Genes usados")
en_n_genes_before=count_genes_in_files("results/predictXcan/test/vcf_1000G_hg37_en", "en_", "en")
en_n_genes_processed=count_genes_in_files_processed("results/predictXcan/test/vcf_1000G_hg37_en/processed", "processed_en_", "proc_en")


n_genes_tissues = pd.merge(mashr_n_genes_before, mashr_n_genes_processed, on='Tejido', how='outer')
n_genes_tissues = pd.merge(n_genes_tissues, en_n_genes_before, on='Tejido', how='outer')
n_genes_tissues = pd.merge(n_genes_tissues, en_n_genes_processed, on='Tejido', how='outer')


n_genes_tissues = n_genes_tissues.sort_values(by='mashr', ascending=False)


# guardamos la tabla en un csv:
n_genes_tissues.to_csv(f"{path}/results/predictXcan/test/n_genes_tissues.csv", index=False)


# contamos numeros de snps tras procesar:
def count_snps(custom_path, file_prefix):
    predict_files = [f for f in os.listdir(custom_path) if f.startswith(file_prefix) and f.endswith("_summary.txt")]
    genes_df = {}

    for predict_file in predict_files:
        start_index_name = predict_file.find(file_prefix) + len(file_prefix)
        end_index_name = predict_file.find("_summary.txt")
        predict_file_name = predict_file[start_index_name:end_index_name]       
        predict_file_path = os.path.join(custom_path, predict_file)
        df_predict = pd.read_csv(predict_file_path, sep="\t")
        n_snps_model=df_predict["n_snps_in_model"].sum()
        n_snps_model_used=df_predict["n_snps_used"].sum()
        genes_df[predict_file_name] = [n_snps_model, n_snps_model_used]
        print(predict_file, "finished")

    genes_df_final = pd.DataFrame.from_dict(genes_df, orient='index', columns=["SNPs en mashr", "SNPs Usados"]).reset_index()
    genes_df_final["SNPs Usados"] = genes_df_final["SNPs Usados"].astype(int)
    genes_df_final = genes_df_final.rename(columns={'index': 'Tejido'})

    return genes_df_final

mashr_snps=count_snps("results/predictXcan/test/vcf_1000G_hg37_mashr", "mashr_")

path_mashr_snps = f"{path}/results/imputado/extracted/all_positions_mashr.txt"

# caculamos el total de SNPs que hay en el modelo cogiendo all_positons_mashr.txt 
with open(path_mashr_snps, 'r') as file:
    line_count = sum(1 for _ in file)
mashr_total_snps = line_count
# el numero de snps usados es el mismo de SNPs que tenemos en nuestro genoma (ya que lo hemos filtrado con los SNPs de mashr), esto lo calculamos 
# en filtrarvcf_position_mashr.sh, es 216665: 

mashr_snps.loc[len(mashr_snps)] = ['Total', mashr_total_snps, 216665]

predixcan_snps_genes = pd.merge(mashr_n_genes_processed, mashr_snps, on='Tejido', how='outer')
predixcan_snps_genes = pd.merge(mashr_n_genes_before, predixcan_snps_genes, on='Tejido', how='outer')

predixcan_snps_genes = predixcan_snps_genes.sort_values(by='Genes en mashr', ascending=False)

predixcan_snps_genes["% Genes usados"]=predixcan_snps_genes["Genes usados"]*100/predixcan_snps_genes["Genes en mashr"]
predixcan_snps_genes["% Genes usados"] = predixcan_snps_genes["% Genes usados"].round(2)

predixcan_snps_genes["% SNPs usados"]=predixcan_snps_genes["SNPs Usados"]*100/predixcan_snps_genes["SNPs en mashr"]
predixcan_snps_genes["% SNPs usados"]= predixcan_snps_genes["% SNPs usados"].round(2)

predixcan_snps_genes = predixcan_snps_genes.sort_values(by='Genes en mashr', ascending=False)

#### ESTO ES LO BUENO
predixcan_snps_genes.to_csv(f"{path}/results/predictXcan/test/n_snps_genes_tissues.csv", index=False) # para cargar la tabla en el informe
# este csv contiene una columna con la misma info que SNPs_tejidos_mashr.csv y SNPs_tejidos.csv 
# (que es redundante ahora, pero habia que comprobar que era igual)


