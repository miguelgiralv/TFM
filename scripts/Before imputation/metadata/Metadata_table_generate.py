import pandas as pd
import GEOparse
import re
import os
import subprocess


############################################
# Generamos una tabla con los metadatos:
def GSM_metadata (gsm_inicial, gsm_final):
    dfs = []
    for i in range(gsm_inicial, gsm_final+1):
        gsm_accession="GSM" + str(i)
        gsm = GEOparse.get_GEO(gsm_accession, destdir= "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GSMs")
        # Extraemos si pasó el control de calidad primero
        characteristics = gsm.metadata["characteristics_ch1"]
        if characteristics[0] == "qc status: Passed QC":
            qc = "Passed"
        elif characteristics[0] == "qc status: Failed QC":
            qc = "Failed"
        else:
            qc = "unknown"    
        # Extraemos el tipo celular de los metadatos
        match = re.search(r'cell type:\s*(\w+)', characteristics[1], re.IGNORECASE)
        if match:
            cell_type = match.group(1)
        else:
            cell_type = "unknown"
        # Extraemos el phenotipo de la enfermedad los metadatos
        if characteristics[2] == "disease status: Unaffected":
            phenotype = "control"
        elif characteristics[2] == "disease status: AD":
            phenotype = "AD"
        else:
            phenotype = "unknown"       
        new_row_df = pd.DataFrame({
            "Individual ID": gsm_accession,
            "QC": qc,
            "Cell type": cell_type,
            "Disease status": phenotype,
            "Individual ID": gsm_accession
        }, index=[0])
        # actualizamos la lista que almacena las filas del df
        dfs.append(new_row_df)
    # generamos el df final
    result_df = pd.concat(dfs, ignore_index=True)
    print(result_df)
    return(result_df)
     
# test
df0=GSM_metadata(894502, 894503)
#generamos varios dfs (algunos gsm ids en medio no son del mismo experimento)
df1=GSM_metadata(894502, 894572).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/metadata_1.txt", index=False, header=True, sep="\t")
df2=GSM_metadata(894589, 894801).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/metadata_2.txt", index=False, header=True, sep="\t")
df3=GSM_metadata(894825, 894895).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/metadata_3.txt", index=False, header=True, sep="\t")
df4=GSM_metadata(894897, 894967).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/metadata_4.txt", index=False, header=True, sep="\t")
df5=GSM_metadata(894980, 895768).to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/metadata_5.txt", index=False, header=True, sep="\t")

# unimos los dfs en uno único que contiene todos los metadatos:
merged_metadata = pd.DataFrame()
metadata_dir = "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/metadatas/"
for metadata_file in os.listdir(metadata_dir):
    file_path = os.path.join(metadata_dir, metadata_file)
    df = pd.read_csv(file_path, delimiter="\t") 
    merged_metadata = pd.concat([merged_metadata, df], ignore_index=True)

#comprobamos la tabla de metadatos:
merged_metadata.head
# vemos si hay columnas redundantes (con un solo valor, tal como cell type=blood) y las eliminamos:
for column_name in merged_metadata.columns:
    if merged_metadata[column_name].nunique() <2:
        merged_metadata=merged_metadata.drop(columns=column_name)
# hemos eliminado cell type porque solo tiene un valor
merged_metadata.head
# guardamos finalmente la tabla con metadatos:
merged_metadata.to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", index=False, header=True, sep="\t")

######################################################################
# Selección de individuos con buena calidad con los metadatos 

# ahora seleccionamos los IDs de los individuos con qc=failed
qc_failed=merged_metadata[merged_metadata["QC"]=="Failed"]
#creamos un df con la lista de ids de muestras a eliminar y su family ID, que es siempre el mismo ()
qc_failed_ids = qc_failed[["Individual ID"]]
qc_failed_ids['Family_id'] = "GSE33528"
# cambiamos el orden para que family_id sea la primera columna
qc_failed_ids = qc_failed_ids[['Family_id', 'Individual ID']]
# lo almacenamos en un txt
qc_failed_ids.to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt", sep="\t", index=False, header=False, lineterminator="\n")

# ahora ejecutamos plink para eliminar las muestras que fallaron el test de calidad y que hemos almacenado en qc_failed_ids.txt
import subprocess
# el comando es:
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528", 
"--remove", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt", 
"--make-bed", "--out", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc"]

# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)

#######################################################################
# Selección de individuos con buena calidad en base al call-rate de variantes (call rate de 0.95) (comando --geno)
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc", 
"--geno", "0.05",
"--out", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr",
"--make-bed"] 
# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)

#######################################################################
# Selección de individuos con buena calidad en base al call-rate de individuos (comando --mind)
plink_command = ["C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe",
"--bfile", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr", 
"--mind", "0.05",
"--out", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr_in",
"--make-bed"] 
# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(plink_command)


#######################################################################
# Control del numero de muestras y variantes que vamos teniendo con el análisis de calidad

#Vemos el número de muestras que teniamos antes del qc:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528.fam", "r") as file:
    row_count_qc_antes = len(file.readlines())    
#y despues del qc
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc.fam", "r") as file:
    row_count_qc_despues = len(file.readlines())    
# y despues de filtrar individuos con call-rate menor:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr.fam", "r") as file:
    row_count_qc_cr_despues = len(file.readlines())   

print("Antes del análisis de calidad inicial, teníamos", row_count_qc_antes, "individuos.","\n",
"Tras el análisis de calidad inicial, nos quedamos con", row_count_qc_despues, "individuos.","\n",
"Hemos eliminado", row_count_qc_antes-row_count_qc_despues, "individuos en el análisis de calidad de los metadatos.","\n",
"Tras el filtrado por call rate, nos quedamos", row_count_qc_cr_despues, "individuos.","\n",
"En este paso de filtrado, hemos eliminado", row_count_qc_despues-row_count_qc_cr_despues, "individuos")

#Vemos el número de variantes que teniamos antes del call-rate:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528.bim", "r") as file:
    row_count_cr_antes = len(file.readlines())    
#y despues de filtrar el call-rate:
with open("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr.bim", "r") as file:
    row_count_cr_despues = len(file.readlines())    
print("Inicialmente teníamos", row_count_cr_antes, "variantes.""\n",
"Tras filtrar variantes con call rate >0.95, nos quedamos", row_count_cr_despues, "variantes","\n",
"Hemos descartado", row_count_cr_antes-row_count_cr_despues, "variantes")

#################################################################################
# mapeo de posiviones desde hg18 a hg38

# mapeamos las coordenadas de hg18 a hg38
# descargamos la cadena de https://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/?C=M;O=A
# Probamos con el liftover de GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard  
#primero generamos el diccionario
liftover_dictionary= ["java", "-jar", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar", "CreateSequenceDictionary",
"REFERENCE=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna",
"OUTPUT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.dict"]
subprocess.run(liftover_dictionary)

#despues ejecutamos liftovervcf
liftover_command = ["java","-Xmx8g","-jar" ,"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar",
"LiftoverVcf",
"I=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_test.vcf",
"O=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/lifted_over.vcf",
"CHAIN=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/hg18ToHg38.over.chain.gz",
"REJECT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/rejected_variants.vcf",
"R=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna"] 

#ejecutamos map_to_UCSC_bed.py:

map_to_bed=["python","C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/map_to_UCSC_bed_python3.py", "-m", 
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map",
"-o", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GSE33528_UCSC"]

subprocess.run(map_to_bed)

# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(liftover_command)

#java -jar -Xmx8g picard.jar LiftoverVcf I=GSE33528_test.vcf O=lifted_over.vcf CHAIN=hg18ToHg38.over.chain.gz REJECT=rejected_variants.vcf R=GCF_000001405.26_GRCh38_genomic.fna

subprocess.run("java -Xmx8g -jar C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar \
C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/LiftoverVcf \
I=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_test.vcf \
O=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/lifted_over.vcf \
CHAIN=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/hg18ToHg38.over.chain.gz \
REJECT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/rejected_variants.vcf \
R=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna", shell=True)

############################################################

# intento mapear de .map a bed usc:

# https://genome.sph.umich.edu/wiki/LiftOver#Method_1
# copio parte del script de https://genome.sph.umich.edu/wiki/LiftMap.py

def myopen(fn):
    import gzip
    try:
        with gzip.open(fn, 'rb') as h:
            ln = h.read(2)  # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn, 'r')
    return gzip.open(fn, 'rt')


def map2bed(fin, fout):
    with open(fout, 'w') as fo:
        for ln in myopen(fin):
            chrom, rs, mdist, pos = ln.split()
            chrom = 'chr' + chrom
            pos = int(pos)
            fo.write('%s\t%d\t%d\t%s\n' % (chrom, pos-1, pos, rs))
    return True
## generamos el archivo bed
map2bed("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_out_ucsc.bed")
 
# ejecutamos liftover en wsl:
#  ./liftOver GSE33528_out_ucsc.bed hg18ToHg38.over.chain.gz output.bed unlifted.bed

