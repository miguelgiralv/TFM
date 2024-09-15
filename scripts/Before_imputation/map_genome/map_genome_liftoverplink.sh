# definimos los paths relativos:
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
work_dir="/home/miguel/liftover"

# Copiaremos los archivos necesarios al directorio de wsl directamente ya que el codigo de liftoverPlink 
# da problemas al tratar con paths de archivos de windows con espacios en medio 
cp "$path/results/plink_data/GSE33528.map" "$work_dir"
cp "$path/software/liftOver" "$work_dir"
cp "$path/software/liftOverPlink.py" "$work_dir"
cp "$path/data/coordinate_map/hg18ToHg19.over.chain.gz" "$work_dir"
mkdir results

#ahora ejecutamos liftoverplink:
python2 "$work_dir/liftOverPlink.py" -m "$work_dir/GSE33528.map" -p "$work_dir/GSE33528.ped" -o "$work_dir/results/output_map" -c "$work_dir/hg18ToHg19.over.chain.gz" -e "$work_dir/liftOver"

# y copiamos los resultados mapeados a su carpeta del repositorio:
cp "$work_dir/results"/* "$path/results/plink_data/classic/processed"

# tras mapear deberemos convertir a plink binario (plink_binary)
