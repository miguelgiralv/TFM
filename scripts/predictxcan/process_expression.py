import os
import pandas as pd

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
### APLICAR A ELASTIC NET

summary_files = [f for f in os.listdir(os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en")) if f.endswith('summary.txt')]

# creamos un diccionario vacio para almacenar los df filtrados 2: Initialize an empty dictionary to store filtered DataFrames
filtered_dataframes = {}

# Step 3: Iterate over each summary file
for summary_file in summary_files:
    # Construct the path to the predict file
    predict_file = summary_file.replace('_summary.txt', '_predict.txt')
    
    # Construct the full path to the summary file and predict file
    summary_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en", summary_file)
    predict_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en", predict_file)
    
    # Step 4: Read the summary file into df_summary
    df_summary = pd.read_csv(summary_file_path, sep="\t")
    
    # Step 5: Identify genes to remove based on 'n_snps_used' column being NA
    genes_to_remove = df_summary.loc[df_summary['n_snps_used'].isna(), 'gene'].tolist()
    
    # Step 6: Read the predict file into df_predict
    df_predict = pd.read_csv(predict_file_path, sep="\t")
    
    # Step 7: Exclude columns from df_predict based on genes_to_remove and 'FID'
    columns_to_keep = [col for col in df_predict.columns if col not in genes_to_remove and col != 'FID']
    
    # Step 8: Filter df_predict to keep only the columns in columns_to_keep
    filtered_df_predict = df_predict[columns_to_keep]
    
    # Step 9: Construct path to save processed file
    path_save = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en/processed", f"processed_{predict_file}")
    
    # Save the filtered DataFrame to CSV
    filtered_df_predict.to_csv(path_save, sep="\t", index=False)
    
    print(f"Processed {predict_file}: Initial dimensions={df_predict.shape}, Filtered dimensions={filtered_df_predict.shape}")
    
    filtered_dataframes[summary_file] = filtered_df_predict

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

### NUMBER OF GENES IN EACH TISSUE FOR ELASTIC NET

predict_files = [f for f in os.listdir(os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en/processed"))]

genes_df_en = {}

for predict_file in predict_files:
    start_index_name = predict_file.find("processed_en_") + len("processed_en_")
    end_index_name = predict_file.find("_predict.txt")
    predict_file_name = predict_file[start_index_name:end_index_name]
      
    predict_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_en/processed", predict_file)
    df_predict = pd.read_csv(predict_file_path, sep="\t")
    n_genes=len(df_predict.columns)-1
    genes_df_en[predict_file_name]=n_genes
    print(predict_file, "finished")

elastic_net_genes = pd.DataFrame(genes_df_en.items(), columns=['Tejido', 'en'])

    
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


predict_files = [f for f in os.listdir(os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr/processed"))]

# creamos un diccionario vacio para almacenar los df filtrados 2: Initialize an empty dictionary to store filtered DataFrames
genes_df_mashr = {}

for predict_file in predict_files:
    start_index_name = predict_file.find("processed_mashr_") + len("processed_mashr_")
    end_index_name = predict_file.find("_predict.txt")
    predict_file_name = predict_file[start_index_name:end_index_name]
      
    predict_file_path = os.path.join(path, "results/predictXcan/test/vcf_1000G_hg37_mashr/processed", predict_file)
    df_predict = pd.read_csv(predict_file_path, sep="\t")
    n_genes=len(df_predict.columns)-1
    genes_df_mashr[predict_file_name]=n_genes
    print(predict_file, "finished")

mashr_genes = pd.DataFrame(genes_df_mashr.items(), columns=['Tejido', 'mashr'])
n_genes_tissues = pd.merge(elastic_net_genes, mashr_genes, on='Tejido', how='outer')
n_genes_tissues.to_csv(f"{path}/results/predictXcan/test/n_genes_tissues.csv", index=False)
