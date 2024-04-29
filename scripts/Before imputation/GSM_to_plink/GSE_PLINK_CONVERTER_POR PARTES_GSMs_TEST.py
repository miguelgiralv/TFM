###########################################
# Definimos funcion para generar los ped

import pandas as pd
import GEOparse
import subprocess

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

###########################################


# convertir los GSM en ped por partes:
# test
GSM_to_plink(894502, 894503).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/GSE33528_0.ped", index=False, header=False, sep=" ")
# el resto:
GSM_to_plink(894502, 894572).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_1.ped", index=False, header=False, sep="\t")
GSM_to_plink(894589, 894719).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_2.ped", index=False, header=False, sep="\t")
GSM_to_plink(894720, 894801).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_3.ped", index=False, header=False, sep="\t")
GSM_to_plink(894825, 894895).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_4.ped", index=False, header=False, sep="\t")
GSM_to_plink(894897, 894920).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_4_2.ped", index=False, header=False, sep="\t")
GSM_to_plink(894921, 894967).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_5.ped", index=False, header=False, sep="\t")
GSM_to_plink(894980, 895122).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_6.ped", index=False, header=False, sep="\t")
GSM_to_plink(895123, 895273).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_7.ped", index=False, header=False, sep="\t")
GSM_to_plink(895274, 895424).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_8.ped", index=False, header=False, sep="\t")
GSM_to_plink(895425, 895575).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_9.ped", index=False, header=False, sep="\t")
GSM_to_plink(895576, 895676).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_10.ped", index=False, header=False, sep="\t")
GSM_to_plink(895677, 895768).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_11.ped", index=False, header=False, sep="\t")

###########################################

# Unimos todos los ped en uno:
def merge_ped_files(input_files, output_file):
    # Abrimos el output file en modo append
    with open(output_file, "a") as outfile:
        for input_file in input_files:
            # recorremos todos los archivos de input para leerlos
            with open(input_file, "r") as infile:
                # escribimos cada linea del inputfile en el output file
                for line in infile:
                    outfile.write(line)
# Creamos la lista de archivos que fusionar
input_files = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_1.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_2.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_3.ped","C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_4.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_4_2.ped", 
               "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_5.ped","C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_6.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_7.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_8.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_9.ped", 
               "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_10.ped", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_11.ped"]

# Creamos la ruta del archivo de output
output_file = "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/merged_ped_file.ped"

# Ejecutamos
merge_ped_files(input_files, output_file)
###########################################
# Cambiamos los "-" en las variantes faltantes por "0" para que no de error convertir a plink binario:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/merged_ped_file.ped", "r") as infile, open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.ped", "w") as outfile:
    # Recorremos todas las linas
    for line in infile:
        # Dividimos las lineas en campo 
        fields = line.split()
        # Iteramos los campos y sustitoumos - por 0
        for i in range(len(fields)):
            if fields[i] == "-":
                fields[i] = "0"
        # Escribimos la linea modificada en el archivo nuevo
        outfile.write("\t".join(fields) + "\n")
        
###########################################
# comprobamos que el numero de muestras coincide con las muestras contando las lineas:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.ped", 'r') as file:
    line_count = sum(1 for line in file)
print(line_count)
# coincide!

# a continuación debemos realizar el mapeo de hg18 a hg38, una vez hecho, generamos el plink binario mapeado:
# Tras esto, ejecutamos en plink en el directorio donde esta el ped y map para generar el plink binario:
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink",
"--file",
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528",
"--make-bed",
"--out",
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/GSE33528"]

subprocess.run(plink_command)


