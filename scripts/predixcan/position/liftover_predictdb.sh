path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# liftover de predictdb 
# hemos creado un archivo que contiene todas las posiciones para liftover (en process_db_mashr_rename.py)
$path/software/liftOver \
    $path/results/imputado/extracted/all_positions_mashr.txt \
    $path/data/predictXcan/hg38ToHg19.over.chain.gz \
    $path/results/imputado/extracted/all_positions_mashr_hg19.txt \
    $path/results/imputado/extracted/all_positions_mashr_hg19_unlifted.txt

# a este archivo de liftover le quitamos chr (chr_remove.py)