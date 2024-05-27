import pandas as pd

######################################################################

path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
# Selección de individuos con buena calidad con los metadatos 
merged_metadata=pd.read_csv(f"{path}/results/metadata_gsm/merged_metadata.txt", sep="\t" )

# ahora seleccionamos los IDs de los individuos con qc=failed
qc_failed=merged_metadata[merged_metadata["QC"]=="Failed"]
#creamos un df con la lista de ids de muestras a eliminar y su family ID, que es siempre el mismo ()
qc_failed_ids = qc_failed[["Individual ID"]]
qc_failed_ids['Family_id'] = "GSE33528"
# cambiamos el orden para que family_id sea la primera columna
qc_failed_ids = qc_failed_ids[['Family_id', 'Individual ID']]
# lo almacenamos en un txt
qc_failed_ids.to_csv(f"{path}/results/metadata_gsm/qc_failed_ids.txt", sep="\t", index=False, header=False, lineterminator="\n")

