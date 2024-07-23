import pandas as pd
import os
import numpy as np


###### buenoo
tissue=predict_files[1]
tissue_pd=pd.DataFrame
path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

files_list = os.listdir(tissue_paths)

predict_files = [file for file in files_list if file.endswith("predict.txt")]
data_frames = []
tissue_names = []
merged_data=pd.DataFrame

for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID", "IID"], axis=1)
    tissue_name = tissue.replace("mashr_", "").replace("_predict.txt", "")
    tissue_pd["tissue"]=tissue_name
    data_frames.append(tissue_pd)
    tissue_names.append(tissue_name)
    print(tissue_name,"cargado")
    
merged_data = pd.concat(data_frames)

print(merged_data.head())

long_data = pd.melt(merged_data.reset_index(), id_vars=['tissue'], var_name='Gene', value_name='Expression')

print(long_data.head())



# Perform ANOVA
genes = merged_data.columns[:-1]  # Exclude the tissue column
anova_results = {}

import scipy.stats as stats
print(merged_data.head())

genes = [genes[0], genes[1]]  # Selects only 2 genes to test

anova_results = {}

for gene in genes:
    # Exclude rows with NA values for the current gene
    non_na_data = merged_data[[gene, 'tissue']].dropna()
    
    # Debugging: Check the data
    print(f"Data for gene {gene}:")
    print(non_na_data.head())
    print(non_na_data['tissue'].value_counts())

    # Ensure that 'tissue' column is a Series
    tissue_series = non_na_data['tissue']
    if isinstance(tissue_series, pd.Series):
        # Get unique tissue types
        tissue_types = tissue_series.unique()
        print(f"Tissue types for gene {gene}: {tissue_types}")

        # Group gene expression levels by tissue type
        expression_levels_by_tissue = [non_na_data[non_na_data['tissue'] == tissue][gene] for tissue in tissue_types]

        for i, levels in enumerate(expression_levels_by_tissue):
            print(f"Tissue {tissue_types[i]} - Number of samples: {len(levels)}")
            print(f"Tissue {tissue_types[i]} - Variance: {np.var(levels)}")
        
        # Perform ANOVA if there are at least two groups to compare
        if len(expression_levels_by_tissue) > 1:
            try:
                anova = stats.f_oneway(*expression_levels_by_tissue)
                anova_results[gene] = anova.pvalue
            except Exception as e:
                print(f"Error for gene {gene}: {e}")
        else:
            anova_results[gene] = None  # No ANOVA performed if less than two groups
    else:
        print(f"'tissue' column is not a Series: {type(tissue_series)}")
    
    print(f"terminado {gene} to go")

# Save results
anova_df = pd.DataFrame(list(anova_results.items()), columns=['Gene', 'p-value'])
anova_df = anova_df.sort_values(by='p-value')
print(anova_df.head())




long=len(genes)
for gene in genes:
    # Exclude rows with NA values for the current gene
    non_na_data = merged_data[[gene, 'tissue']].dropna()
    
    # Group gene expression levels by tissue type
    expression_levels_by_tissue = [non_na_data[non_na_data['tissue'] == tissue][gene] for tissue in non_na_data['tissue'].unique()]
    
    # Perform ANOVA if there are at least two groups to compare
    if len(expression_levels_by_tissue) > 1:
        anova = stats.f_oneway(*expression_levels_by_tissue)
        anova_results[gene] = anova.pvalue
    long=long-1
    print("terminado",gene,long,"to go")

# Save the p-value results in a DataFrame and sort by p-value
anova_df = pd.DataFrame(list(anova_results.items()), columns=['Gene', 'p-value'])
anova_df = anova_df.sort_values(by='p-value')

print(anova_df.head())


# Optional: Save the DataFrame to a CSV file
anova_df.to_csv('anova_pvalues.csv', index=False)

# Print the ordered DataFrame
print(anova_df)









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
path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

files_list = os.listdir(tissue_paths)

predict_files = [file for file in files_list if file.endswith("predict.txt")]
data_frames = []
tissue_names = []
# quitar

for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID", "IID"], axis=1)
    means = tissue_pd.mean()
    means_df = pd.DataFrame([means])
    tissue_name = tissue.replace("mashr_", "").replace("_predict.txt", "")
    means_df["tissue"]=tissue_name
    data_frames.append(means_df)
    tissue_names.append(tissue_name)
    print(tissue_name,"cargado")
    
merged_data = pd.concat(data_frames)
print(merged_data.head())

print(merged_data["tissue"])

genes = merged_data.columns[:-1]
anova_results = {}

tissue=tissue_names[1]
gene=genes[1]

for gene in genes:
    merged_data[gene]
    non_na_data = merged_data[[gene, 'tissue']].dropna()



merged_data = pd.concat(data_frames)
print(merged_data.head())
print(merged_data["tissue"])

# Perform ANOVA
genes = merged_data.columns[:-1]  # Exclude the tissue column
anova_results = {}

# Save the p-value results in a DataFrame and sort by p-value
anova_df = pd.DataFrame(list(anova_results.items()), columns=['Gene', 'p-value'])
anova_df = anova_df.sort_values(by='p-value')

# Optional: Save the DataFrame to a CSV file
anova_df.to_csv('anova_pvalues.csv', index=False)

# Print the ordered DataFrame
print(anova_df)



    for tissue in tissue_names:
        
        merged_data[merged_data["tissue"]== tissue]
    anova = stats.f_oneway(*[merged_data[merged_data['tissue'] == tissue][gene] for tissue in merged_data['tissue'].unique()])
    anova_results[gene] = anova.pvalue

for gene in genes:
    # Initialize an empty list to hold the data for each tissue
    groups = []
    # Group by tissue and get the values for the current gene
    for tissue, group in merged_data.groupby('tissue'):
        groups.append(group[gene].values)
    # Perform ANOVA
    anova = stats.f_oneway(*groups)
    # Store the p-value in the dictionary
    anova_results[gene] = anova.pvalue
    
    
# Step 4: Save the p-value results in a DataFrame
anova_df = pd.DataFrame(list(anova_results.items()), columns=['Gene', 'p-value'])

# Step 5: Order the DataFrame from smallest to largest p-value
anova_df = anova_df.sort_values(by='p-value')

# Optional: Save the DataFrame to a CSV file
anova_df.to_csv('anova_pvalues.csv', index=False)

# Print the ordered DataFrame
print(anova_df)



tissue=predict_files[1]
tissue_pd=pd.DataFrame
path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

files_list = os.listdir(tissue_paths)

predict_files = [file for file in files_list if file.endswith("predict.txt")]
data_frames = []

for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID", "IID"], axis=1)
    tissue_pd=tissue_pd.T
    # quitamos el indice que esta ahora como headers
    tissue_pd.columns = tissue_pd.iloc[0]
    tissue_pd = tissue_pd[1:]
    means = tissue_pd.mean()
    means_df = pd.DataFrame([means])
    tissue_name = tissue.replace("mashr_", "").replace("predict.txt", "")
    means_df["tissue"]=tissue_name
    data_frames.append(means_df)
    # Renombramos las columnas para que cada gen una tenga el nombre del tejido
    tissue_pd.columns = [f"{col}_{tissue_name}" if col != 'IID' else col for col in tissue_pd.columns]
    data_frames.append(means_df)
    
    print(tissue,"cargado")
  
merged_data = pd.concat(data_frames, axis=1, join='inner')
print(merged_data.head())

import scipy.stats as stats

def clean_column_name(column_name):
    """
    Clean and standardize the column name.
    Extracts the part of the column name after the second underscore,
    converts it to lowercase, and removes any trailing underscores.
    """
    parts = column_name.split('_')
    if len(parts) > 2:
        tissue_name = '_'.join(parts[2:]).rstrip('_')
        tissue_name = tissue_name.lower()
        tissue_name = tissue_name.strip()
        return tissue_name
    return column_name

# Clean all column names
cleaned_columns = [clean_column_name(col) for col in merged_data.columns]

# Store cleaned column names as tissue_names
tissue_names = set(cleaned_columns)

def perform_anova(expression_values):
    """Perform one-way ANOVA."""
    if len(expression_values) < 2:
        return None, None  # Need at least two groups
    try:
        f_stat, p_value = stats.f_oneway(*expression_values)
        return f_stat, p_value
    except Exception as e:
        print(f"ANOVA error: {e}")
        return None, None

# Prepare DataFrame to store ANOVA results
anova_results = pd.DataFrame(columns=['Gene', 'F-statistic', 'p-value'])

# List to hold results before concatenation
results_list = []

# Perform ANOVA for each gene
for gene in merged_data.index:
    gene_data = merged_data.loc[gene]
    # Extract expression values for each tissue
    tissue_values = []
    for tissue in tissue_names:
        cols = [col for col in merged_data.columns if tissue in col]
        values = gene_data[cols].values
        if len(values) > 1:  # Ensure there are enough data points
            tissue_values.append(values)
        else:
            print(f"Not enough data for tissue {tissue} for gene {gene}")
    if len(tissue_values) > 1:  # Ensure there are at least two groups to compare
        f_stat, p_value = perform_anova(tissue_values)
        if f_stat is not None:
            results_list.append({'Gene': gene, 'F-statistic': f_stat, 'p-value': p_value})
    else:
        print(f"Not enough groups to perform ANOVA for gene {gene}")

# Convert results list to DataFrame and concatenate
anova_results = pd.concat([pd.DataFrame(results_list)], ignore_index=True)

# Print the result
print(anova_results)