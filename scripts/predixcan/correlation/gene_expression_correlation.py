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
# emtemos un contador para tener todos los genes distintos
number = 1
#  Creamos un diccionario vacio para almacenar una matriz de correlacion para cada gen
all_correlations = {}  
# Iteramos por todos los dfs de expresion genica entodos  los tejidos
for df_index, df in enumerate(data_frames):
    # Loop por los genes/columnas del tejido
    for col_name in df.columns:
        if col_name != 'IID' and col_name != 'tissue' and col_name not in all_correlations:
            # df vacio para m,eter un gen
            gene_df = pd.DataFrame()
            # cogemos la expresion de ese gen en todos los tejidos
            for i, tissue_df in enumerate(data_frames):
                tissue_name = tissue_names[i]
                if col_name in tissue_df.columns:
                    gene_tissue = f"{col_name}_{tissue_name}"
                    gene_df[gene_tissue] = tissue_df[col_name]            
            # Matriz de correlacion para el df
            if not gene_df.empty:
                correlation_matrix = gene_df.corr()
                new_index = [name.split('_', 1)[-1] for name in correlation_matrix.index]
                correlation_matrix.index = new_index
                if col_name not in all_correlations:
                    all_correlations[col_name] = []
                # lo metemos en el  diccionario con el nomnbre del gen
                all_correlations[col_name].append(correlation_matrix)
        
            number += 1
            print(f"Processed gene {col_name}, count: {number}")


flattened_correlations = {}

for gene, matrices in all_correlations.items():
    if len(matrices) > 1:
        flattened_correlations[gene] = pd.concat(matrices, axis=0)
    else:
        flattened_correlations[gene] = matrices[0]

merged_df = pd.concat(flattened_correlations.values(), axis=1)

transposed_df = merged_df.T

# COMPROBAR QUE ES LO QUE QUERIA CARLOS Y LUEGO HACER LOOP POR TODOS LOS TEJIDOS PARA NO DEJARME NINGUN GEN, LUEGO HACER MERGE QUITANDO LOS GENES REDUNDANTES
# pendiente
merged_df = pd.concat(all_correlations, axis=1)
transposed_df = merged_df.T
transposed_df.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation.txt", sep="\t", header= True)

transposed_df=pd.read_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation.txt", sep="\t")

def split_gene_tissue(value):
    parts = value.split('_', 1)
    gene = parts[0]
    tissue = parts[1] if len(parts) > 1 else None
    return gene, tissue

transposed_df[['Gene', 'Tissue']] = transposed_df['Unnamed: 0'].apply(lambda x: pd.Series(split_gene_tissue(x)))

transposed_df.drop(columns=['Unnamed: 0'], inplace=True)

tissues = transposed_df.columns[1:-2]  # Adjusting for the 'Gene' and 'Tissue' columns

mean_correlation_matrix = pd.DataFrame(index=tissues, columns=tissues)

for tissue1 in tissues:
    for tissue2 in tissues:
        tissue1_data = transposed_df[transposed_df['Tissue'] == tissue1]
        # Extraemos las correlaciones 
        correlations = tissue1_data[tissue2]
        # Calculamos correlacion media entre tissue 1 y 2
        mean_correlation_matrix.loc[tissue1, tissue2] = correlations.mean()

mean_correlation_matrix = mean_correlation_matrix.astype(float)

mean_correlation_matrix.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_mean_matrix.txt", sep="\t", header= True)

mean_correlation_matrix=pd.read_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_mean_matrix.txt", sep="\t")

mean_correlation_matrix.set_index('Unnamed: 0', inplace=True)

sorted_columns = sorted(mean_correlation_matrix.columns)
mean_correlation_matrix = mean_correlation_matrix[sorted_columns]
mean_correlation_matrix = mean_correlation_matrix.sort_index()

plt.figure(figsize=(30, 30))
ax=sns.heatmap(mean_correlation_matrix, annot=False, cmap='Reds', vmin=0, vmax=1,square=True)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xticks(rotation=90)  
plt.yticks(rotation=0)  
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
# creamos una funcion para quitar todos los genes que se expresan en ese tejido para quitar los valores con correlacion 1 pq es el mismo gen
def replace_values(row, columns):
    for col in columns:
        if row['Unnamed: 0'].endswith(f'_{col}'):
            row[col] = pd.NA
    return row

# Aplicamos a cada fila
transposed_df2
transposed_df2 = transposed_df2.apply(lambda row: replace_values(row, columns), axis=1)
transposed_df2_filter = transposed_df2.drop(columns=['Unnamed: 0'])
transposed_df2_filter = transposed_df2_filter[sorted(transposed_df2_filter.columns)]
transposed_df2_melt = transposed_df2_filter.melt(var_name='Tissue', value_name='Correlation')

transposed_df2_melt.to_csv(f"{path}/results/predictXcan/test/vcf_1000G_hg37_mashr/correlation_boxplot.txt", sep="\t", header= True)

# y lo representamos en un boxplot
plt.figure(figsize=(12, 8))
sns.boxplot(x='Tissue', y='Correlation', data=transposed_df2_melt)
plt.xticks(rotation=90)  
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
plt.xticks(rotation=90)  
plt.yticks(rotation=0)   
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

# comprobamos que las columnas existen
existing_columns = [col for col in selected_columns if col in mean_correlation_matrix.columns]
print("Columns selected:", existing_columns)

# Seleccionamos las columnas de los tejidos selected
mean_correlation_matrix_filtered = mean_correlation_matrix[['Unnamed: 0'] + existing_columns]

filtered_df = mean_correlation_matrix_filtered[mean_correlation_matrix_filtered['Unnamed: 0'].isin(selected_columns)]

filtered_df.set_index('Unnamed: 0', inplace=True)
filtered_df = filtered_df.sort_index()

plt.figure(figsize=(12, 10))
ax=sns.heatmap(filtered_df, annot=False, cmap='Reds', vmin=0, vmax=1)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
plt.xticks(rotation=90)  
plt.yticks(rotation=0)   
plt.xlabel("Tissues", fontsize=18, labelpad=30)
plt.ylabel("Tissues", fontsize=18, labelpad=30)
plt.show()
