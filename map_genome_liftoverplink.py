
###########################################
# Definimos funcion para generar los ped

import pandas as pd
import GEOparse

# Columnas de los SNP en la tabla de cada muestra
selected_columns = ["Allele1-Top", "Allele2-Top"]
# Lista para meter a los objetos de df
# Cargamos el archivo soft de varios GSM
def GSM_to_plink (gsm_inicial, gsm_final):
    dfs = []
    for i in range(gsm_inicial, gsm_final+1):
        gsm_accession="GSM" + str(i)
        gsm_test = GEOparse.get_GEO(gsm_accession, destdir= "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GSMs")
        # ELIMINART ESTA PARTE, VALE CON gsm_test_id
        # Extraemos el phenotipo de los metadatos
        characteristics = gsm_test.metadata["characteristics_ch1"]
        if characteristics[2] == "disease status: Unaffected":
            phenotype_value = 1
        elif characteristics[2] == "disease status: AD":
            phenotype_value = 2
        else:
            phenotype_value = -9
        # Extraemos el sexo comparando un par de sondas del cromosoma Y:
        Y_variant_not_present = gsm_test.table[(gsm_test.table["ID_REF"].isin(["rs2032590-124_T_R_IFB1141653252:", "rs2032597-124_B_R_IFB1141671058:"]) & (gsm_test.table["Allele1-Top"] == '-'))]
        if len(Y_variant_not_present)>1:
            # si la tabla de variantes de Y tiene más de 1 elemento, significa que las dos variantes de Y 
            # están y tienen como valor "-", por lo tanto es una mujer (sex = 2 en ped)
            sexo = 2
        else:
            # si la tabla no tiene tiene 2 variantes asociadas a Y con valor "-",será un hombre (sex = 1 en ped)
            sexo = 1
        # Construimos un df que solo contiene las columnas de SNPs de la gsm.table
        selected_df = gsm_test.table[selected_columns]
        # Creamos un df vacio de una sola fila, con los datos de la muestra:
        new_row_df = pd.DataFrame({
            "Family ID": "GSE33528",
            "Male parent": 0,
            "Female parent": 0,
            "Sex": sexo,
            "Phenotype": phenotype_value,
            "Individual ID": gsm_accession
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

    return(result_df)

###########################################

# Generamos el archivo map
gpl = GEOparse.get_GEO('GPL14932',
                           destdir="C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO")
selected_columns = ["Chr","Name", "MapInfo"]
map_table = gpl.table[selected_columns]
map_table["position"] = 0
order_columns= ["Chr","Name","position","MapInfo"]
map_table.loc[:,order_columns].to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/GSE33528_2.map", sep='\t', index=False, header=False)
map_table["Name"].to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/GSE33528_rsID.txt", sep='\t', index=False, header=False)




