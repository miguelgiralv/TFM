
# Copiaremos los archivos necesarios al directorio de wsl directamente ya que el codigo de liftoverPlink da problemas al tratar con paths de archivos con espacios en medio 
cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map" /home/miguel/liftover

cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/liftOver" /home/miguel/liftover

cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/liftOverPlink.py" /home/miguel/liftover

cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.ped" /home/miguel/liftover

cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/hg18ToHg38.over.chain.gz" /home/miguel/liftover

#ahora ejecutamos liftoverplink:
python2 liftOverPlink.py -m /home/miguel/liftover/GSE33528.map -p /home/miguel/liftover/GSE33528.ped -o /home/miguel/liftover/results/output_map -c /home/miguel/liftover/hg18ToHg38.over.chain.gz -e /home/miguel/liftover/liftOver

# y copiamos los resultados mapeados a su carpeta:
cp /home/miguel/liftover/results/* "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/processed"


