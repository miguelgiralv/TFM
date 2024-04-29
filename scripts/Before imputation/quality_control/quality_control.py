import pandas as pd
import subprocess

######################################################################
# Selección de individuos con buena calidad con los metadatos 

merged_metadata=pd.read_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/merged_metadata.txt", sep="\t" )

# ahora seleccionamos los IDs de los individuos con qc=failed
qc_failed=merged_metadata[merged_metadata["QC"]=="Failed"]
#creamos un df con la lista de ids de muestras a eliminar y su family ID, que es siempre el mismo ()
qc_failed_ids = qc_failed[["Individual ID"]]
qc_failed_ids['Family_id'] = "GSE33528"
# cambiamos el orden para que family_id sea la primera columna
qc_failed_ids = qc_failed_ids[['Family_id', 'Individual ID']]
# lo almacenamos en un txt
qc_failed_ids.to_csv("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt", sep="\t", index=False, header=False, lineterminator="\n")

# ahora ejecutaremos plink para eliminar las muestras que fallaron el test de calidad y que hemos almacenado en qc_failed_ids.txt (plink_QC.bat)
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

