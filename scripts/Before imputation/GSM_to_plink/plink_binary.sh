
# Una vez hecho el liftover a hg19, generamos el plink binario mapeado con los archivos map y ped mapeados (output_map):
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM"

"$path/software/plink" --file "$path/results/plink_data/classic/processed/output_map" \
--allow-extra-chr --make-bed --out "$path/results/plink_data/binary/raw/GSE33528"

# Finalmente corregiremos el sexo de nuestro plink binario con la funcion check-sex:
"$path/software/plink" --bfile "$path/results/plink_data/binary/raw/GSE33528" \
--check-sex --allow-extra-chr --out "$path/results/metadata_gsm/GSE33528_sex_check_results"

# Aplicamos la predicción de plink a nuestro archivo plink binario con la funcion impute-sex para tener el sexo corregido:
"$path/software/plink" --bfile "$path/results/plink_data/binary/raw/GSE33528" \
--impute-sex --allow-extra-chr --make-bed --out \
"$path/results/plink_data/binary/processed/GSE33528_sex"
