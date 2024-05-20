path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM"

# Ejecutaremos plink para eliminar las muestras que fallaron el test de calidad según los metadatos y que hemos almacenado en qc_failed_ids.txt (en el 
# script de python). El comando es:
"$path/software/plink" --bfile "$path/results/plink_data/binary/processed/GSE33528_sex" --remove "$path/results/metadata_gsm/qc_failed_ids.txt" \
--make-bed --out "$path/results/plink_data/binary/processed/GSE33528_sex_qc"

#######################################################################
# Seleccionamos ahora las variantes con buena calidad en base al call-rate de variantes (call rate de 0.95) (comando --geno)
"$path/software/plink" --bfile "$path/results/plink_data/binary/processed/GSE33528_sex_qc" --geno 0.05 \
--out "$path/results/plink_data/binary/processed/GSE33528_qc_cr" --make-bed

#######################################################################
# Finalmente, hacemos una selección de individuos con buena calidad en base al call-rate de individuos (comando --mind)
"$path/software/plink" --bfile "$path/results/plink_data/binary/processed/GSE33528_qc_cr" --mind 0.05 --out "$path/results/plink_data/binary/processed/GSE33528_qc_cr_in" --make-bed 



# el pipeline continua en report samples
