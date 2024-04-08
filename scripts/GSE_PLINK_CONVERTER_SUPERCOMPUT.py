import pandas as pd
from GEOparse import GEOparse
import os

# Columnas de los SNP en la tabla de cada muestra
selected_columns = ["Allele1-Top", "Allele2-Top"]
# Lista para meter a los objetos de df
dfs = []
# Cargamos el archivo soft de GSE
gse_test = GEOparse.get_GEO(filepath="C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/samples/GSE33528_family.soft.gz")
# Iteramos por todas las muestras
def geo_gsm_to_plink_partes(number_1, number_2):
    for gsm_test_id, gsm_test in gse_test.gsms.items(): 
        # Convertimos el identificador individual de la muestra en número
        numerical_part = gsm_test_id.replace("GSM", "")
        gsm_number = int(numerical_part)
        # De esta forma podemos seleccionar las muestras que vamos a convertir en .ped y podemos convertir 
        # todas las muestras por partes
        if number_1 <= gsm_number <= number_2:
            # ELIMINART ESTA PARTE, VALE CON gsm_test_id
            accession = gsm_test.metadata["geo_accession"]
            # Extraemos el phenotipo de los metadatos
            characteristics = gsm_test.metadata["characteristics_ch1"]
            if characteristics[2] == "disease status: Unaffected":
                phenotype_value = 1
            elif characteristics[2] == "disease status: Affected":
                phenotype_value = 2
            else:
                phenotype_value = -9
            # Construimos un df que solo contiene las columnas de SNPs de la gsm.table
            selected_df = gsm_test.table[selected_columns]

            # Creamos un df vacio de una sola fila, con los datos de la muestra:
            new_row_df = pd.DataFrame({
                "Family ID": "GSE33528",
                "Male parent": 0,
                "Female parent": 0,
                "Sex": "other",
                "Phenotype": phenotype_value,
                "Individual ID": accession[0]
            }, index=[0])

            # Concatenamos las columnas de SNPs (alelo1 y alelo2) de la muestra en una sola fila para las columnas del .ped
            score_theta_values = selected_df.values.flatten()
            score_theta_df = pd.DataFrame(score_theta_values).T
            # Concatenamos la nueva fila con info de la muestra con el df que contiene los snps en una sola fila
            new_row_df = pd.concat([new_row_df, score_theta_df], axis=1)
            # Añadimos la nueva fila del df en la lista vacia que se ha creado
            dfs.append(new_row_df)

    # Concatenaremos todas las filas del df que estaban almacenadas en la lista dfs
    result_df = pd.concat(dfs, ignore_index=True)

    orden_columnas = ["Family ID", "Individual ID", "Male parent", "Female parent", "Sex", "Phenotype"]
    columnas_snp_orden = sorted(result_df.columns.difference(orden_columnas))
    final_columns_order = orden_columnas + list(columnas_snp_orden)
    result_df = result_df[final_columns_order]

    print(result_df)
    directory = "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/"
    output_file = os.path.join(directory, f"GSE33528_{number_1}_{number_2}.ped")
    result_df.to_csv(output_file, index=False)

