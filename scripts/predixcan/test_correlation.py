import pandas as pd
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns



path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

files_list = os.listdir(tissue_paths)

predict_files = [file for file in files_list if file.endswith("predict.txt")]
data_frames = []
tissue_names = []
merged_data=pd.DataFrame

for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID"], axis=1)
    tissue_name = tissue.replace("mashr_", "").replace("_predict.txt", "")
    data_frames.append(tissue_pd)
    tissue_names.append(tissue_name)
    print(tissue_name,"cargado")

data_frames[1]
###bien

genes=[]
for tissue in len(data_frames):
    tissue_df=data_frames[tissue]
    for col in tissue_df.columns:
        if col != 'IID':
            col=[]
            for tissue in len(data_frames):
                col={}
                
                
###BIEN! COMPROBAR CON CARLOS Y HACER EL LOOP POR TODOS LOS TEJIDOS PARA NO DEJARME NINGUN GEN, LUEGO HACER MERGE EN UNA                
number=1                
all_correlations = []
# Iterate through columns of the first DataFrame (excluding 'IID')
for col_name in data_frames[1].columns:
    if col_name != 'IID' and col_name != 'tissue':
        # Create a DataFrame to hold data for the current gene across all tissues
        gene_df = pd.DataFrame()
        for i, df in enumerate(data_frames):
            tissue_name = tissue_names[i]
            if col_name in df.columns:
                gene_tissue = f"{col_name}_{tissue_name}"
                gene_df[gene_tissue] = df[col_name]       
        # Calculate the correlation matrix for the current gene
        if not gene_df.empty:
            correlation_matrix = gene_df.corr()
            new_index = [name.split('_', 1)[-1] for name in correlation_matrix.index]
            correlation_matrix.index = new_index
        all_correlations.append(correlation_matrix)
        number=number+1
        print(number)

df = df.loc[:, ~df.columns.duplicated()]

# COMPROBAR QUE ES LO QUE QUERIA CARLOS Y LUEGO HACER LOOP POR TODOS LOS TEJIDOS PARA NO DEJARME NINGUN GEN, LUEGO HACER MERGE QUITANDO LOS GENES REDUNDANTES
merged_df = pd.concat(all_correlations, axis=1)

transposed_df = merged_df.T
transposed_df.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation.txt", sep="\t", header= True)

transposed_df=pd.read_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation.txt", sep="\t")

def split_gene_tissue(value):
    # Split only at the first underscore
    parts = value.split('_', 1)
    gene = parts[0]
    tissue = parts[1] if len(parts) > 1 else None
    return gene, tissue

# Apply the function to split the column into 'Gene' and 'Tissue'
transposed_df[['Gene', 'Tissue']] = transposed_df['Unnamed: 0'].apply(lambda x: pd.Series(split_gene_tissue(x)))

transposed_df.drop(columns=['Unnamed: 0'], inplace=True)

tissues = transposed_df.columns[1:-2]  # Adjusting for the 'Gene' and 'Tissue' columns

mean_correlation_matrix = pd.DataFrame(index=tissues, columns=tissues)

for tissue1 in tissues:
    for tissue2 in tissues:
        # Filter the data for the first tissue
        tissue1_data = transposed_df[transposed_df['Tissue'] == tissue1]
        
        # Extract the correlations for tissue2
        correlations = tissue1_data[tissue2]
        
        # Calculate the mean correlation for the gene in tissue1 with its expression in tissue2
        mean_correlation_matrix.loc[tissue1, tissue2] = correlations.mean()

mean_correlation_matrix = mean_correlation_matrix.astype(float)

mean_correlation_matrix.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_mean_matrix.txt", sep="\t", header= True)

mean_correlation_matrix=pd.read_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_mean_matrix.txt", sep="\t")

mean_correlation_matrix.set_index('Unnamed: 0', inplace=True)
# Reorder columns alphabetically
sorted_columns = sorted(mean_correlation_matrix.columns)
mean_correlation_matrix = mean_correlation_matrix[sorted_columns]
# Sort rows by index (if needed)
mean_correlation_matrix = mean_correlation_matrix.sort_index()

plt.figure(figsize=(30, 30))
ax=sns.heatmap(mean_correlation_matrix, annot=False, cmap='Reds', vmin=0, vmax=1,square=True)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xticks(rotation=90)  # Rotate x-axis labels if needed for better readability
plt.yticks(rotation=0)   # Keep y-axis labels horizontal
plt.xlabel("Tissues", fontsize=18, labelpad=40)
plt.ylabel("Tissues", fontsize=18, labelpad=40)
plt.show()


#solo con tejido de cerebro
brain_tissues = [tissue for tissue in mean_correlation_matrix.columns if tissue.lower().startswith('brain')]

filtered_matrix = mean_correlation_matrix.loc[brain_tissues, brain_tissues]

plt.figure(figsize=(12, 10))
ax=sns.heatmap(filtered_matrix, annot=False, cmap='Reds', vmin=0, vmax=1)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xticks(rotation=90)  # Rotate x-axis labels if needed for better readability
plt.yticks(rotation=0)   # Keep y-axis labels horizontal
plt.xlabel("Tissues", fontsize=18, labelpad=30)
plt.ylabel("Tissues", fontsize=18, labelpad=30)
plt.show()

# solo con 19 tejidos:
mean_correlation_matrix=pd.read_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_mean_matrix.txt", sep="\t")

selected_columns = [
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Substantia_nigra",
    "Cells_EBV-transformed_lymphocytes",
    "Heart_Atrial_Appendage",
    "Heart_Left_Ventricle",
    "Whole_Blood"
]

# Check if the columns exist in the DataFrame
existing_columns = [col for col in selected_columns if col in mean_correlation_matrix.columns]
print("Columns selected:", existing_columns)

# Filter the DataFrame to include only the selected columns
mean_correlation_matrix_filtered = mean_correlation_matrix[['Unnamed: 0'] + existing_columns]

filtered_df = mean_correlation_matrix_filtered[mean_correlation_matrix_filtered['Unnamed: 0'].isin(selected_columns)]

filtered_df.set_index('Unnamed: 0', inplace=True)
filtered_df = filtered_df.sort_index()

# Ensure the 'Unnamed: 0' column is set as the index for the heatmap
plt.figure(figsize=(12, 10))
ax=sns.heatmap(filtered_df, annot=False, cmap='Reds', vmin=0, vmax=1)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xticks(rotation=90)  # Rotate x-axis labels if needed for better readability
plt.yticks(rotation=0)   # Keep y-axis labels horizontal
plt.xlabel("Tissues", fontsize=18, labelpad=30)
plt.ylabel("Tissues", fontsize=18, labelpad=30)
plt.show()


### no valido elimianr

    
for tissue in predict_files:
    tissue_pd=pd.read_csv(os.path.join(tissue_paths, tissue), sep="\t",  on_bad_lines='skip')
    tissue_pd=tissue_pd.drop(["FID"], axis=1)
    tissue_name = tissue.replace("mashr_", "").replace("_predict.txt", "")
    new_columns = {}
    for col in tissue_pd.columns:
        if col != 'IID':  
            new_columns[col] = f"{col}_{tissue_name}"
        else:
            new_columns[col] = col  
    tissue_pd.rename(columns=new_columns, inplace=True)  
    data_frames.append(tissue_pd)
    tissue_names.append(tissue_name)
    print(tissue_name,"cargado")



merged_data = pd.concat(data_frames)




tissue_paths = f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

# List files in directory
files_list = os.listdir(tissue_paths)

# Filter files that end with "predict.txt"
predict_files = [file for file in files_list if file.endswith("predict.txt")]

# Initialize a variable to store the merged DataFrame progressively
merged_data = None

# Process and merge each file progressively
for tissue in predict_files:
    file_path = os.path.join(tissue_paths, tissue)
    
    try:
        # Read the file into a DataFrame
        tissue_pd = pd.read_csv(file_path, sep="\t", on_bad_lines='skip')
        
        # Drop the 'FID' column if it exists
        tissue_pd = tissue_pd.drop(["FID"], axis=1, errors='ignore')
        
        # Extract tissue name from file name
        tissue_name = tissue.replace("mashr_", "").replace("_predict.txt", "")
        
        # Rename columns
        new_columns = {col: f"{col}_{tissue_name}" if col != 'IID' else col for col in tissue_pd.columns}
        tissue_pd.rename(columns=new_columns, inplace=True)
        
        # Print status
        print(f"{tissue_name} loaded")
        
        # Merge the current DataFrame with the accumulated DataFrame
        if merged_data is None:
            merged_data = tissue_pd
        else:
            merged_data = pd.concat([merged_data, tissue_pd], ignore_index=True)       
    except Exception as e:
        print(f"Error processing {tissue}: {e}")

merged_data = pd.concat(data_frames)

# hacer contraste columna a columna:
for col in merged_data.columns:
    col
    



print(merged_data.head())

# todas las columnas menos el tejido
genes = merged_data.columns[:-1]  
genes = merged_data.columns[0:3]
anova_results = {}

long = len(genes)

for gene in genes:
    non_na_data = merged_data[[gene, 'tissue']].dropna()

    # Group gene expression levels by tissue type
    expression_levels_by_tissue = [non_na_data[non_na_data['tissue'] == tissue][gene] for tissue in non_na_data['tissue'].unique()]
    
    # Perform ANOVA if there are at least two groups to compare
    if len(expression_levels_by_tissue) > 1:
        try:
            anova = stats.f_oneway(*expression_levels_by_tissue)
            anova_results[gene] = anova.pvalue
        except Exception as e:
            print(f"Error for gene {gene}: {e}")
    else:
        anova_results[gene] = None  # No ANOVA performed if less than two groups
    
    long -= 1
    print(f"terminado {gene}, {long} to go")

# Save results
anova_df = pd.DataFrame(list(anova_results.items()), columns=['Gene', 'p-value'])
anova_df = anova_df.sort_values(by='p-value')
print(anova_df.head())
