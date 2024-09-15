import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

tissue_paths=f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr"

# cargamos los archivos de expresion de predixcan (predict.txt)
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

# comprobamos que se ha cargado bien
data_frames[1]
                        
###BIEN! COMPROBAR CON CARLOS Y HACER EL LOOP POR TODOS LOS TEJIDOS PARA NO DEJARME NINGUN GEN, LUEGO HACER MERGE EN UNA                
number=1                
all_correlations = []
# Iteramos por todas las columnas de los genes 
for element in len(data_frames):
    for col_name in data_frames[1].columns:
        if col_name != 'IID' and col_name != 'tissue':
            # iteramos por todos los genes y almacenamos los datos de cada gen en todos los tejidos en un df
            gene_df = pd.DataFrame()
            for i, df in enumerate(data_frames):
                tissue_name = tissue_names[i]
                if col_name in df.columns:
                    gene_tissue = f"{col_name}_{tissue_name}"
                    gene_df[gene_tissue] = df[col_name]       
            # calculamos entonces la matriz de correlacion para ese gen
            if not gene_df.empty:
                correlation_matrix = gene_df.corr()
                new_index = [name.split('_', 1)[-1] for name in correlation_matrix.index]
                correlation_matrix.index = new_index
            all_correlations.append(correlation_matrix)
            number=number+1
            print(number)


df = df.loc[:, ~df.columns.duplicated()]

############################ ALTERNATIVA CON TODOS LOS
number = 1
all_correlations = {}  # Initialize as a dictionary to store correlation matrices for each gene
# Iterate through each data frame in data_frames
for df_index, df in enumerate(data_frames):
    # Iterate through the columns (genes) of the current data frame
    for col_name in df.columns:
        if col_name != 'IID' and col_name != 'tissue' and col_name not in all_correlations:
            # Create an empty dataframe to store gene data across tissues
            gene_df = pd.DataFrame()
            # Iterate over all data frames and collect data for this gene
            for i, tissue_df in enumerate(data_frames):
                tissue_name = tissue_names[i]
                if col_name in tissue_df.columns:
                    gene_tissue = f"{col_name}_{tissue_name}"
                    gene_df[gene_tissue] = tissue_df[col_name]            
            # Calculate the correlation matrix for this gene
            if not gene_df.empty:
                correlation_matrix = gene_df.corr()
                # Rename the index to just tissue names
                new_index = [name.split('_', 1)[-1] for name in correlation_matrix.index]
                correlation_matrix.index = new_index
                # Initialize the key for the gene if it doesn't exist in all_correlations
                if col_name not in all_correlations:
                    all_correlations[col_name] = []
                # Append the correlation matrix for this gene
                all_correlations[col_name].append(correlation_matrix)
            # Increment and print progress
            number += 1
            print(f"Processed gene {col_name}, count: {number}")


























number = 1
all_correlations = []
for df_index, df in enumerate(data_frames):
    # Iterate over all columns of genes in the current DataFrame
    for col_name in df.columns:
        if col_name != 'IID' and col_name != 'tissue':
            # Create a DataFrame to store data for the current gene across all tissues
            gene_df = pd.DataFrame()            
            for i, tissue_df in enumerate(data_frames):
                tissue_name = tissue_names[i]
                if col_name in tissue_df.columns:
                    gene_tissue = f"{col_name}_{tissue_name}"
                    gene_df[gene_tissue] = tissue_df[col_name]            
            # Calculate the correlation matrix for the current gene
            if not gene_df.empty:
                correlation_matrix = gene_df.corr()
                new_index = [name.split('_', 1)[-1] for name in correlation_matrix.index]
                correlation_matrix.index = new_index
                correlation_matrix.columns = new_index
                all_correlations.append(correlation_matrix)
            number += 1
            print(number)

# Combine all correlation matrices into a single DataFrame
# Use MultiIndex for columns to keep track of different genes
combined_correlations = pd.concat(all_correlations, axis=1, keys=[f"Gene_{i}" for i in range(len(all_correlations))])

len(combined_correlations)

# Drop duplicate columns and rows
# Flatten the MultiIndex columns to a single level if needed
combined_correlations.columns = combined_correlations.columns.map('_'.join)
combined_correlations = combined_correlations.loc[~combined_correlations.index.duplicated(keep='first')]
combined_correlations = combined_correlations.loc[:, ~combined_correlations.columns.duplicated(keep='first')]

print("Combined correlations matrix with duplicates dropped:")
print(combined_correlations)

combined_correlations.head()



# COMPROBAR QUE ES LO QUE QUERIA CARLOS Y LUEGO HACER LOOP POR TODOS LOS TEJIDOS PARA NO DEJARME NINGUN GEN, LUEGO HACER MERGE QUITANDO LOS GENES REDUNDANTES
# pendiente
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



##### calculamos la media de correlacion entre tejidos:
mean_correlation_matrix_filter = mean_correlation_matrix.replace(1, pd.NA)
mean_correlation_matrix_filter = mean_correlation_matrix_filter.drop(columns=['Unnamed: 0'])
mean_correlation_matrix_melt = mean_correlation_matrix_filter.melt(var_name='Tissue', value_name='Correlation')
print(mean_correlation_matrix_melt.head())

# eliminamos los genes que se correlacionan con su mismo tejido (da 1)

transposed_df2=transposed_df
columns = transposed_df2.columns[transposed_df2.columns != 'Unnamed: 0']
# Function to replace values
def replace_values(row, columns):
    for col in columns:
        if row['Unnamed: 0'].endswith(f'_{col}'):
            row[col] = pd.NA
    return row
# Apply the function to each row
transposed_df2 = transposed_df2.apply(lambda row: replace_values(row, columns), axis=1)
transposed_df2_filter = transposed_df2.drop(columns=['Unnamed: 0'])
transposed_df2_filter = transposed_df2_filter[sorted(transposed_df2_filter.columns)]
transposed_df2_melt = transposed_df2_filter.melt(var_name='Tissue', value_name='Correlation')

transposed_df2_melt.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_boxplot.txt", sep="\t", header= True)

# y lo representamos en un boxplot
plt.figure(figsize=(12, 8))
sns.boxplot(x='Tissue', y='Correlation', data=transposed_df2_melt)
plt.xticks(rotation=90)  # Rotate x labels for better readability
plt.title('Correlation of gene expression accross Tissues', fontsize = 20)
plt.xlabel("Tissues", fontsize=16)
plt.ylabel("Correlation", fontsize=18)
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
